#include "stereo-calibration.h"

#include <exception>


stereo_calibration::
stereo_calibration(const std::string& _filename)
{
  read_xml(_filename);
}


stereo_calibration::
stereo_calibration(
    const std::vector<cv::Mat>& _images_1,
    const std::vector<cv::Mat>& _images_2,
    const size_t _x,
    const size_t _y,
    const double _square_size
)
  : stereo_calibration(
      camera_calibration(_images_1, _x, _y, _square_size),
      camera_calibration(_images_2, _x, _y, _square_size),
      _images_1,
      _images_2,
      _x, _y, _square_size)
{ }


stereo_calibration::
stereo_calibration(
    const camera_calibration& _camera_1,
    const camera_calibration& _camera_2,
    const std::vector<cv::Mat>& _images_1,
    const std::vector<cv::Mat>& _images_2,
    const size_t _x,
    const size_t _y,
    const double _square_size
)
  : camera_1(_camera_1), camera_2(_camera_2)
{
  // Assert that we got the same number of images from each camera.
  if(_images_1.size() != _images_2.size())
    throw std::runtime_error("stereo_calibration error: need same number of "
        "images from each camera.");

  std::cout << "Computing stereo calibration...\n";

  // To calibrate, we need one point in both image and world coords for each
  // intersection of four squares.
  const size_t num_x = _x - 1,
               num_y = _y - 1;

  // Detect the square intersections in each image.
  std::vector<std::vector<cv::Vec2f>> image_1_points, image_2_points;
  for(size_t i = 0; i < _images_1.size(); ++i) {
    try {
      // Take both corresponding images or neither of them.
      auto points_1 = extract_image_points(_images_1[i], num_x, num_y);
      auto points_2 = extract_image_points(_images_2[i], num_x, num_y);
      image_1_points.push_back(std::move(points_1));
      image_2_points.push_back(std::move(points_2));
    }
    catch(std::exception& _e) {
      std::cout << "Error in processing image pair " << i << ":\n\t" << _e.what()
                << "\nSkipping this image pair.\n"
                << std::endl;
    }
  }

  // Assert that we got at least two usable images.
  if(image_1_points.size() < 2) {
    std::cout << "\nAt least two usable images are required to calibrate."
              << std::endl;
    std::exit(0);
  }

  // Generate world points for the target.
  auto target_points = create_world_points(num_x, num_y, _square_size);

  // We need one set of these for each calibration image.
  std::vector<std::vector<cv::Vec3f>> world_points;
  for(size_t i = 0; i < image_1_points.size(); ++i)
    world_points.push_back(target_points);

  // Set output image size.
  image_size = camera_1.image_size;

  // Compute the stereo calibration.
#if CV_VERSION_MAJOR == 3
  const double error = cv::stereoCalibrate(
      world_points, image_1_points, image_2_points,
      camera_1.intrinsics, camera_1.distortion,
      camera_2.intrinsics, camera_2.distortion,
      image_size,
      R, T, E, F,
      cv::CALIB_FIX_INTRINSIC,
      cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 300,
        1e-8)
  );
#else
  const double error = cv::stereoCalibrate(
      world_points, image_1_points, image_2_points,
      camera_1.intrinsics, camera_1.distortion,
      camera_2.intrinsics, camera_2.distortion,
      image_size,
      R, T, E, F,
      cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 300,
        1e-8),
      cv::CALIB_FIX_INTRINSIC
  );
#endif

  std::cout << "Error: " << error << std::endl;
}


void
stereo_calibration::
compute_rectify_maps()
{
  // Make space for the rectification parameters.
  cv::Mat R1, R2, Q;

  // Rectify the calibration.
  cv::stereoRectify(
      camera_1.intrinsics, camera_1.distortion,
      camera_2.intrinsics, camera_2.distortion,
      image_size, R, T,
      R1, R2, P1, P2, Q,
      CV_CALIB_ZERO_DISPARITY,
      0,
      image_size,
      &rect_1, &rect_2);

  // Compute undistort maps.
  maps.resize(4);
  cv::initUndistortRectifyMap(camera_1.intrinsics, camera_1.distortion, R1, P1,
      image_size, CV_32FC1, maps[0], maps[1]);

  cv::initUndistortRectifyMap(camera_2.intrinsics, camera_2.distortion, R2, P2,
      image_size, CV_32FC1, maps[2], maps[3]);
}

/*----------------------------------------------------------------------------*/

void
stereo_calibration::
read_xml(
    const std::string& _filename
) {
  cv::FileStorage fs(_filename, cv::FileStorage::READ);

  // Read the sizes into a buffer to avoid problems with older OpenCV (v 2.4.5).
  std::vector<int> buffer;

  fs["intrinsics-1"] >> camera_1.intrinsics;
  fs["distortion-1"] >> camera_1.distortion;
  //fs["image_size-1"] >> camera_1.image_size;
  fs["image_size-1"] >> buffer;
  camera_1.image_size = cv::Size(buffer[0], buffer[1]);

  fs["intrinsics-2"] >> camera_2.intrinsics;
  fs["distortion-2"] >> camera_2.distortion;
  //fs["image_size-2"] >> camera_2.image_size;
  buffer.clear();
  fs["image_size-2"] >> buffer;
  camera_2.image_size = cv::Size(buffer[0], buffer[1]);

  fs["E"] >> E;
  fs["F"] >> F;
  fs["R"] >> R;
  fs["T"] >> T;

  fs["P1"] >> P1;
  fs["P2"] >> P2;

  buffer.clear();
  fs["image_size"] >> buffer;
  image_size = cv::Size(buffer[0], buffer[1]);

  // Same problems with rectangles and old OpenCV.
  //fs["rect_1"] >> rect_1;
  //fs["rect_2"] >> rect_2;
  buffer.clear();
  fs["rect_1"] >> buffer;
  rect_1 = cv::Rect(buffer[0], buffer[1], buffer[2], buffer[3]);
  buffer.clear();
  fs["rect_2"] >> buffer;
  rect_2 = cv::Rect(buffer[0], buffer[1], buffer[2], buffer[3]);

  fs.release();
}


void
stereo_calibration::
write_xml(
    const std::string& _filename
) {
  cv::FileStorage fs(_filename, cv::FileStorage::WRITE);

  fs << "intrinsics-1" << camera_1.intrinsics
     << "distortion-1" << camera_1.distortion
     << "image_size-1" << camera_1.image_size;

  fs << "intrinsics-2" << camera_2.intrinsics
     << "distortion-2" << camera_2.distortion
     << "image_size-2" << camera_2.image_size;

  fs << "E" << E
     << "F" << F
     << "R" << R
     << "T" << T;

  fs << "P1" << P1
     << "P2" << P2;

  fs << "image_size" << image_size
     << "rect_1" << rect_1
     << "rect_2" << rect_2;

  fs.release();
}
