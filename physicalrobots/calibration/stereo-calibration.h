#ifndef _STEREO_CALIBRATION_H_
#define _STEREO_CALIBRATION_H_

#include <string>
#include <vector>

#include "calibration.h"

#include <opencv2/opencv.hpp>


////////////////////////////////////////////////////////////////////////////////
/// Calibration for a stereo camera pair with fixed inter-camera transformation.
////////////////////////////////////////////////////////////////////////////////
struct stereo_calibration
{

  ///@name Internal State
  ///@{

  /// Individual calibrations for each camera.
  camera_calibration camera_1, camera_2;

  /// The essential matrix, fundamental matrix, and relative
  /// rotation/translation.
  cv::Mat E, F, R, T;

  /// The stereo-rectified projection matrices.
  cv::Mat P1, P2;

  /// Undistortion maps for each camera.
  std::vector<cv::Mat> maps;

  /// Size that is usable in both images.
  cv::Size image_size;

  /// Rectangles of usable pixels within each image.
  cv::Rect rect_1, rect_2;

  ///@}
  ///@name Construction
  ///@{

  /// Create a stereo calibration from a saved file.
  /// @param _filename The name of the saved XML file.
  stereo_calibration(const std::string& _filename);

  /// Calibrate a stereo rig, starting from existing calibrations for the
  /// individual cameras. Then determine the stereo relationships.
  /// @param _camera_1 The existing calibration for the first camera.
  /// @param _camera_2 The existing calibration for the second camera.
  /// @param _images_1 The images from the first camera.
  /// @param _images_2 The images from the first camera.
  /// @param _x The number of squares in the target width-wise.
  /// @param _y The number of squares in the target height-wise.
  /// @param _square_size The size of the target squares.
  stereo_calibration(
      const camera_calibration& _camera_1,
      const camera_calibration& _camera_2,
      const std::vector<cv::Mat>& _images_1,
      const std::vector<cv::Mat>& _images_2,
      const size_t _x,
      const size_t _y,
      const double _square_size);

  /// Calibrate a stereo rig all in one go.
  /// @param _camera_1 The existing calibration for the first camera.
  /// @param _camera_2 The existing calibration for the second camera.
  /// @param _images_1 The images from the first camera.
  /// @param _images_2 The images from the first camera.
  /// @param _x The number of squares in the target width-wise.
  /// @param _y The number of squares in the target height-wise.
  /// @param _square_size The size of the target squares.
  stereo_calibration(
      const std::vector<cv::Mat>& _images_1,
      const std::vector<cv::Mat>& _images_2,
      const size_t _x,
      const size_t _y,
      const double _square_size);

  /// Compute the rectification for the camera pair.
  void compute_rectify_maps();

  ///@}
  ///@name IO
  ///@{

  /// Read a calibration from an XML file.
  /// @param _filename The XML file name.
  void read_xml(const std::string& _filename);

  /// Write this calibration to an XML file.
  /// @param _filename The XML file name.
  void write_xml(const std::string& _filename);

  ///@}

};

#endif
