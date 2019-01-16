#ifndef MARKER_PERCEPT_H_
#define MARKER_PERCEPT_H_

#include <iostream>
#include <vector>

namespace aruco {
  class Marker;
}


////////////////////////////////////////////////////////////////////////////////
/// Describes the perception of an ArUco marker, giving id, distance, and
/// heading from camera.
///
/// The heading measured here is the angle theta w.r.t. the robot's forward
/// direction:
///
///        forward
///           ^
///     +theta|
///        \--|
///         \ |
///          \|
///  <--------|-------> right
///
////////////////////////////////////////////////////////////////////////////////
class MarkerPercept {

  ///@name Internal State
  ///@{

  int m_id;          ///< The marker's ID.
  double m_distance; ///< The marker's distance from the camera center.
  double m_heading;  ///< The marker's bearing relative to camera view direction.
  double m_angle;    ///< The marker's angle relative to -camera view direction.

  ///@}

  public:

    ///@name Construction
    ///@{

    MarkerPercept(const aruco::Marker& _m);

    ///@}
    ///@name Accessors
    ///@{

    int ID() const noexcept;
    double Distance() const noexcept;
    double Heading() const noexcept;
    double Angle() const noexcept;

    /// Compute the camera's position in the marker's local frame (where +X is
    /// the surface normal, +Z is up, and +Y is stage-right when looking at the
    /// marker).
    std::vector<double> GetCameraPositionInMarkerFrame() const;

    ///@}
    ///@name Equality
    ///@{

    bool operator==(const MarkerPercept& _a) const noexcept;
    bool operator!=(const MarkerPercept& _a) const noexcept;

    ///@}

  private:

    ///@name Measurement Helpers
    ///@{

    /// Compute the distance from the camera to a given marker.
    /// @param _m The marker to localize.
    /// @return The distance to the marker.
    double ComputeDistance(const aruco::Marker& _m) const;

    /// Compute the angular heading from the camera to a given marker.
    /// @param _m The marker to localize.
    /// @return The angular heading to the marker.
    double ComputeHeading(const aruco::Marker& _m) const;

    /// Compute the marker angle as seen in the camera.
    /// The value will be in the range (-pi / 2, pi / 2), with angle increasing
    /// as the marker is rotated in the up direction.
    /// @param _m The marker to localize.
    /// @return The angle of the marker.
    double ComputeAngle(const aruco::Marker& _m) const;

    ///@}

};

/// Output stream overload for printing percepts to cout/cerr.
std::ostream& operator<<(std::ostream& _os, const MarkerPercept& _p);

#endif
