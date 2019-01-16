#ifndef ROBOT_CONTROLLER_H_
#define ROBOT_CONTROLLER_H_

#include <atomic>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>


namespace PlayerCc
{
  class PlayerClient;
  class Position2dProxy;
  class BumperProxy;
}

namespace aruco
{
  class Marker;
};

const double PI = 3.14;


////////////////////////////////////////////////////////////////////////////////
/// A simple controller for sending commands to an iCreate robot.
////////////////////////////////////////////////////////////////////////////////
class RobotController {

  private:

    ///@name Local Types
    ///@{

    /// A queued command for the create to execute.
    struct Command
    {
      const double translation; ///< The translation velocity.
      const double rotation;    ///< The rotation velocity.
      const double time;        ///< The execution length.
    };

    ///@}
    ///@name Internal State
    ///@{

    // Player components
    PlayerCc::PlayerClient* m_client{nullptr};        ///< Player client object.
    PlayerCc::Position2dProxy* m_position2d{nullptr}; ///< Position control.
    PlayerCc::BumperProxy* m_bumper{nullptr};         ///< Bump sensor.

    // Command queue
    std::thread m_commandThread;        ///< Thread for managing command queue.
    std::mutex m_commandQueueLock;      ///< Lock for the command queue.
    std::queue<Command> m_commandQueue; ///< Queue of commands to be executed.
    volatile std::atomic<bool> m_commandThreadRunning{false};

    ///@}

  public:

    ///@name Construction
    ///@{

    /// Construct an iCreate controller.
    /// @param _ip The IP address of the create's netbook.
    RobotController(const std::string& _ip);

    ~RobotController();

    ///@}
    ///@name Command Queue
    ///@{

    /// Add a new command to the robot's queue.
    /// @param _translation The translational velocity.
    /// @param _rotation The rotational velocity.
    /// @param _time The time to execute the command in seconds.
    void EnqueueCommand(const double _translation, const double _rotation,
        const double _time);

    /// Start the command thread and send commands to the robot.
    void StartCommandQueue();

    ///@}
    ///@name Movement Functions
    ///@{

    /// MoveToPoint, which uses Rotate and Translate, attempts to taxi the robot
    /// from its current coordinates to new ones on a straight line.
    ///
    /// However, because the Create has very imprecise tracking of its rotation,
    /// this problem becomes difficult. There is a threshold turn rate (~0.15)
    /// below: which no turns will register in the odometry. Because the Create
    /// calculates its odometry based on absolute distance traveled and amount
    /// turned, its perceived x, y, theta can get really bad, really fast.
    ///
    /// The strategy that MoveToPoint attempts is to always stay above the lower
    /// threshold for turn rate (0.15) and just modulate the duration of the
    /// turn. This is still imprecise, but less so. The robot will also ignore a
    /// heading error unless it becomes some delta away from where it needs to
    /// be (hardcoded at the moment). When its heading error becomes delta, the
    /// robot will attempt to turn again and continue toward its destination.
    ///
    /// Because the odometry of the create is so bad, any real world application
    /// that requires even semi-precise point-point movement requires the robot
    /// be able to localize (ideally using landmarks like markers). Maybe in the
    /// future, one can add a digital compass sensor to the Create, and then
    /// manage odometry on the client side more accurately.
    void MoveToPoint(const double _x, const double _y,
        const double _epsilon = .1);

    void Rotate(double _rads);

    void Translate(double _meters);

    ///@}
    ///@name Odometry Functions
    ///@{

    std::vector<double> GetOdometry();

    void SetOdometry(const double _x, const double _y, const double _yaw);

    ///@}

  private:

    /// Set the speed of the robot through the client and sleep the thread to
    /// allow it to move at that speed for a set ammount of time.
    /// @param _linear The linear speed.
    /// @param _angular The angular speed.
    /// @param _microseconds The time to execute the command, in microseconds.
    /// @note After the command ends, the robot will not zero its speed. Call
    ///       this function with zero'd parameters to stop the robot.
    void SetSpeedAndUSleep(const double _linear, const double _angular,
        const unsigned int _microseconds);

    /// Pull fresh sensor data from the robot.
    void Read();

};

// Revisit later when we are ready to integrate input into the simulator.
void GetXYFromMarker(aruco::Marker& _markers);

#endif
