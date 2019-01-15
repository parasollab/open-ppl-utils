#include "MarkerPercept.h"

#include "aruco.h"

#include <ctgmath>
#include <iomanip>

#ifndef PI
#define PI 3.1415926
#endif


/*---------------------------- Construction ----------------------------------*/

MarkerPercept::
MarkerPercept(const aruco::Marker& _m) {
  m_id = _m.id;
  m_distance = ComputeDistance(_m);
  m_heading  = ComputeHeading(_m);
  m_angle    = ComputeAngle(_m);

  // For debugging.
#if 0
  std::cout << "Marker " << m_id << " data:"
            << "\nx: " << _m.Tvec.ptr<float>(0)[0]
            << "\ny: " << _m.Tvec.ptr<float>(0)[1]
            << "\nz: " << _m.Tvec.ptr<float>(0)[2]
            << "\nrx: " << _m.Rvec.ptr<float>(0)[0]
            << "\nry: " << _m.Rvec.ptr<float>(0)[1]
            << "\nrz: " << _m.Rvec.ptr<float>(0)[2]
            << std::endl;
#endif
}

/*-------------------------------- Accessors ---------------------------------*/

int
MarkerPercept::
ID() const noexcept {
  return m_id;
}


double
MarkerPercept::
Distance() const noexcept {
  return m_distance;
}


double
MarkerPercept::
Heading() const noexcept {
  return m_heading;
}


double
MarkerPercept::
Angle() const noexcept {
  return m_angle;
}


std::vector<double>
MarkerPercept::
GetCameraPositionInMarkerFrame() const {
  // Determine the camera's X, Y, and theta coordinates. Theta is the rotation
  // of the camera's viewing direction from the marker's +X axis, where a
  // positive rotation is about +Z (same for camera and marker frames).

  const double phi = PI/2 + m_angle;
  const double psi = phi + m_heading;

  const double x = m_distance * std::cos(psi),
               y = m_distance * std::sin(psi),
               t = phi + PI;

  return {x, y, t};
}

/*------------------------------ Equality ------------------------------------*/

bool
MarkerPercept::
operator==(const MarkerPercept& _p) const noexcept {
  return _p.m_id == m_id
     and _p.m_distance == m_distance
     and _p.m_heading == m_heading;
}


bool
MarkerPercept::
operator!=(const MarkerPercept& _p) const noexcept {
  return !(*this == _p);
}

/*-------------------------- Measurement Helpers -----------------------------*/

double
MarkerPercept::
ComputeDistance(const aruco::Marker& _m) const {
  // Do not use the Y distance because we assume markers are on the plane.
  return std::sqrt(std::pow(_m.Tvec.ptr<float>(0)[0], 2)
                 //+ std::pow(_m.Tvec.ptr<float>(0)[1], 2)
                 + std::pow(_m.Tvec.ptr<float>(0)[2], 2));
}


double
MarkerPercept::
ComputeHeading(const aruco::Marker& _m) const {
  // Aruco +Z is forward, +X is to the right, +Y is down into the plane.
  const double forwardDistance = _m.Tvec.ptr<float>(0)[2],
               rightDistance   = _m.Tvec.ptr<float>(0)[0];

  // Bearing should be positive if marker is to left of principal
  // axis and negative if marker is to right of principal axis (0 for marker on
  // principle axis).
  return std::atan2(forwardDistance, rightDistance) - PI / 2;
}


double
MarkerPercept::
ComputeAngle(const aruco::Marker& _m) const {
  // The angle frame is strange. When situated at the marker center and looking
  // outward in the facing direction, +X is right, +Y is out, +Z is up.
  // However, aruco uses [-pi:pi] instead of [0:2pi], so the x rotation will
  // change sign when the camera moves above/below the marker.
  const bool flipSign = _m.Rvec.ptr<float>(0)[0] < 0;
  return -_m.Rvec.ptr<float>(0)[2] * (flipSign ? -1 : 1);
}

/*----------------------------------------------------------------------------*/

std::ostream&
operator<<(std::ostream& _os, const MarkerPercept& _p) {
  _os << "{id: " << _p.ID()
      << ", distance: " << std::setprecision(3) << _p.Distance()
      << ", heading: " << std::setprecision(3) << _p.Heading()
      << "}";
  return _os;
}
