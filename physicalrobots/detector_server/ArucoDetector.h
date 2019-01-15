#ifndef ARUCO_DETECTOR_H_
#define ARUCO_DETECTOR_H_

#include "MarkerPercept.h"

#include "aruco.h"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include "cvdrawingutils.h"
#include "cvdrawingutils.h"

#include <string>
#include <vector>


////////////////////////////////////////////////////////////////////////////////
/// Manages detection of ArUco markers from a webcam. Assumes that the robot
/// moves along the ground plane, and that markers are oriented upright with
/// their face normals parallel to the plane.
////////////////////////////////////////////////////////////////////////////////
class ArucoDetector {

  public:

    ///@name Construction
    ///@{

    /// Construct a detector object.
    /// @param _filename The calibration filename.
    ArucoDetector(const std::string& _filename);

    ///@}
    ///@name Interface
    ///@{

    /// Try to detect any markers in the input stream.
    /// @return A vector of marker percepts describing the observation.
    std::vector<MarkerPercept> operator()();

    /// Toggle display of the input video stream.
    /// @param _d True enables display, false disables.
    void Display(const bool _d = true);

    ///@}

  private:

    ///@name Capture Helpers
    ///@{

    /// Capture and store the current image from the video stream.
    /// @param _skipFrames The number of buffered frames to skip before
    ///                    capturing the next image. We think the buffer holds
    ///                    about 5-8 frames. This should be 0 if you are
    ///                    continuously streaming from the camera.
    void CaptureVideo(const size_t _skipFrames = 8);

    /// Get the visible markers from the last captured image.
    /// @return A vector of marker data.
    std::vector<aruco::Marker> GetVisibleMarkers();

    ///@}
    ///@name Display Helpers
    ///@{

    /// Show the captured and thresholded images. The markers are shown in the
    /// captured image.
    /// @param _markers The markers to display.
    void UpdateDisplay(const std::vector<aruco::Marker>& _markers);

    /// Draw a set of markers into an image.
    /// @param _markers The markers to draw.
    /// @param _image The image to draw them on.
    void DrawMarkers(const std::vector<aruco::Marker>& _markers, cv::Mat& _image);

    ///@}
    ///@name Setup Helpers
    ///@{

    /// Read a calibration XML file produced with our calibration tool.
    /// @param _filename The calibration file name.
    void ReadCalibrationFile(const std::string& _filename);

    ///@}
    ///@name Internal State
    ///@

    cv::VideoCapture m_capture;               ///< The video capture stream.
    cv::Mat m_image;                          ///< Last captured image.
    aruco::CameraParameters m_cameraParams;   ///< The camera properties.
    aruco::MarkerDetector m_detector;         ///< Marker detection object.
    bool m_display{false};                    ///< Show video display?

    ///@}

};

#endif
