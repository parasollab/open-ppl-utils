#include "RobotController.h"

#include <opencv2/calib3d/calib3d.hpp>
#include "../aruco/src/aruco.h"
#include "libplayerc++/playerc++.h"

#include <chrono>

#include <unistd.h>

/*------------------------------ Construction --------------------------------*/

RobotController::
RobotController(const std::string& _ip)
  : m_client(new PlayerCc::PlayerClient(_ip, 6665)),
    m_position2d(new PlayerCc::Position2dProxy(m_client, 0))
{
  m_client->StartThread();
  Read();
}


RobotController::
~RobotController()
{
  m_position2d->SetSpeed(0, 0);

  m_client->Stop();

  m_commandThreadRunning = false;
  if(m_commandThread.joinable())
    m_commandThread.join();

  delete m_position2d;
  delete m_client;
}

/*------------------------------ Command Queue -------------------------------*/

void
RobotController::
EnqueueCommand(const double _translation, const double _rotation,
    const double _time)
{
  m_commandQueueLock.lock();
  m_commandQueue.push({_translation, _rotation, _time});
  m_commandQueueLock.unlock();
}


void
RobotController::
StartCommandQueue()
{
  // Guard against trying to start the queue again while it's already running.
  if(m_commandThreadRunning)
    return;
  m_commandThreadRunning = true;

  // Create a functor to execute queued commands. While the queue is empty, this
  // will wait for 100 ms before checking again.
  auto workFunction = [this]()
  {
    while(this->m_commandThreadRunning)
    {
      m_commandQueueLock.lock();

      // If there is no command to execute, create a 'wait' command for 100 ms.
      if(this->m_commandQueue.empty())
        this->m_commandQueue.push({0, 0, .1});

      // Get next command from the queue.
      Command next = this->m_commandQueue.front();
      this->m_commandQueue.pop();

      m_commandQueueLock.unlock();

      this->m_position2d->SetSpeed(next.translation, next.rotation);
      usleep(next.time * 1e6);
    }
  };

  // Start the queue thread.
  m_commandThread = std::thread(workFunction);
}

/*--------------------------- Movement Functions -----------------------------*/

void
RobotController::
MoveToPoint(const double _x, const double _y, const double _epsilon)
{
  while(true)
  {
    auto odom = GetOdometry();
    const double dx = (_x - odom[0]);
    const double dy = (_y - odom[1]);
    const double linearDistance = std::sqrt(dx*dx + dy*dy);

    // The robot is close enough to the point when it's within _epsilon units.
    if(linearDistance <= _epsilon)
      return;

    // Compute angular difference normalized to [-PI, PI].
    double angularDistance = atan2(dy, dx) - odom[2];
    if(angularDistance > PI)
      angularDistance -= 2 * PI;
    else if(angularDistance < -PI)
      angularDistance += 2 * PI;

    // If we need to rotate more than 10 degrees, rotate, else translate.
    if(std::fabs(angularDistance) >= DTOR(10))
      Rotate(angularDistance);
    else
      Translate(std::min(linearDistance, .15));
  }
}


void
RobotController::
Rotate(double _rads)
{
  SetSpeedAndUSleep(0.35, _rads, 0);
  return;
}


void
RobotController::
Translate(double _meters)
{
  // We can only move so accurately, and this is probably a very generous
  // estimation of how much.
  if(fabs(_meters) < 0.02)
    return;

  double sleepTime = 25000.0;
  double moveSpeed = 0.25;
  if(_meters < 0.0)
    moveSpeed *= -1.0;

  int microsecondsToTranslate = (int)(fabs (_meters*1e6 / moveSpeed));
  for(int secondsLeft = microsecondsToTranslate;
      secondsLeft > 0; secondsLeft -= sleepTime) {
    SetSpeedAndUSleep(moveSpeed, 0.0, sleepTime);
  }
  SetSpeedAndUSleep(0.0,0.0,0);
  return;
}

/*--------------------------- Odometry Functions -----------------------------*/

std::vector<double>
RobotController::
GetOdometry()
{
  return {m_position2d->GetXPos(),
          m_position2d->GetYPos(),
          m_position2d->GetYaw()};
}


void
RobotController::
SetOdometry(const double _x, const double _y, const double _yaw)
{
  m_position2d->SetOdometry(_x, _y, _yaw);
}

/*-------------------------------- Helpers -----------------------------------*/

void
RobotController::
SetSpeedAndUSleep(const double _x, const double _a,
    const unsigned int _microseconds)
{
  m_position2d->SetSpeed(_x, _a/2.00);
  usleep(_microseconds);
}


void
RobotController::
Read()
{
  m_client->Read();
}

/*----------------------------------------------------------------------------*/


void
GetXYFromMarker(aruco::Marker& _markers)
{
  cv::Mat rVec = _markers.Rvec;
  // This output shows the amount of rotation about the axis of rotation based
  // on the rodriguez vector
  // double norm = sqrt(rVec.at<float>(0)*rVec.at<float>(0) +
  //  rVec.at<float>(1)*rVec.at<float>(1) +
  // rVec.at<float>(2)*rVec.at<float>(2));
  // cout << "RodriguesMag: " << RTOD(norm) << endl;

  //convert to a roation matrix
  cv::Mat rotMat;
  cv::Rodrigues(rVec, rotMat);

  // theta is rotation about the markers y-axis
  // phi is rotation about the markers x-axis
  // psi is gold and is roation about the markers z-axis
  const double theta = -1 * std::asin(rotMat.at<float>(0, 2));
  //const double phi = asin(rotMat.at<float>(1, 2)/ std::cos(theta));
  const double psi = acos(rotMat.at<float>(0, 0) / std::cos(theta));

  //marker is not considered rotated if it is facing towards -x, the
  //marker's positive z is "up" and positive "x" is left
  //r is distance to the marker
  const double tempX = _markers.Tvec.at<float>(2,0) ;
  const double tempY = _markers.Tvec.at<float>(0,0) * -1.0;
  const double r = sqrt(pow(tempX,2) + pow(tempY,2));

  //alpha will be rotation for the camera to align with the marker
  const double alpha = atan2(tempY, tempX);

  //gamma is global rotation of marker
  double gamma = 0;

  //rHeading is the robots global heading
  double sign = _markers.Rvec.at<float>(1,0) < 0 ? 1.0 : -1.0;
  double rHeading = gamma + (PI-psi) * (sign);

  // Find the heading from the robot's position towards the marker.
  double mHeading = rHeading + alpha;

  const double x = 1 - r * std::cos(mHeading);
  const double y = -r * std::sin(mHeading);

  std::cout << "Current position " << x << " " <<  y << " " << rHeading << endl;
}
