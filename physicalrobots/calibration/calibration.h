#ifndef _CALIBRATION_H_
#define _CALIBRATION_H_

#include <cstddef>
#include <vector>

#include <opencv2/opencv.hpp>


/// Extract image points from a calibration image.
/// @param _image The calibration image.
/// @param _w The number of inner corners in the target width-wise.
/// @param _h The number of inner corners in the target height-wise.
/// @return A set of extracted points.
std::vector<cv::Vec2f>
extract_image_points(
    const cv::Mat& _image,
    const size_t _w,
    const size_t _h
);


/// Create world points for a calibration image.
/// @param _w The number of inner corners width-wise.
/// @param _h The number of inner corners height-wise.
/// @param _square_size The length of a square side.
/// @return A set of world points for the target assuming the upper-left inner
///         corner will be (0,0,0), the entire target lies in the plane z = 0,
///         and +x aims right and +y aims down.
std::vector<cv::Vec3f>
create_world_points(
    const size_t _w,
    const size_t _h,
    const double _square_size
);


////////////////////////////////////////////////////////////////////////////////
/// A complete camera calibration with IO support.
////////////////////////////////////////////////////////////////////////////////
struct camera_calibration
{

  ///@name Internal State
  ///@{

  cv::Mat_<float> intrinsics, distortion;

  cv::Size image_size;

  std::vector<cv::Mat> remaps;

  ///@}
  ///@name Construction
  ///@{

  /// Create an empty calibration.
  camera_calibration();

  /// Create a calibration from intrinsics and distortion parameters.
  /// @param _intrinsics The camera intrinsics.
  /// @param _distortion The distortion parameters.
  /// @param _image_size The size of the images taken by this camera.
  camera_calibration(const cv::Mat& _intrinsics, const cv::Mat& _distortion,
      const cv::Size& _image_size);

  /// Load a calibration from an XML file.
  /// @param _filename The XML file name.
  camera_calibration(const std::string& _filename);

  /// Automatically calibrate a camera using images of a checkerboard target.
  /// @param _images The target images.
  /// @param _x The number of squares across the width of the target.
  /// @param _y The number of squares across the height of the target.
  /// @param _square_size The length of a checker square side.
  camera_calibration(
      const std::vector<cv::Mat>& _images,
      const size_t _x,
      const size_t _y,
      const double _square_size
  );


  void compute_rectify_maps();

  ///@}
  ///@name I/O
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
