#include "calibration.h"

#include <exception>
#include <iostream>


std::vector<cv::Vec2f>
extract_image_points(
    const cv::Mat& _image,
    const size_t _w,
    const size_t _h
) {
  // Create a size object to describe the number of internal corners.
  const cv::Size pattern_size(_w, _h);

  // Convert image to grey scale.
  cv::Mat grey_scale;
  cvtColor(_image, grey_scale, CV_BGR2GRAY);

  // Extract corner positions.
  cv::Mat corners;
  const bool success = cv::findChessboardCorners(grey_scale, pattern_size,
      corners);

  // Assert success.
  if(!success)
    throw std::runtime_error("findChessboardCorners failed.");

  // Sub-pixel refinement.
  const cv::Size refinement_window(5,5), dead_zone(1,1);
  cv::cornerSubPix(grey_scale, corners, refinement_window, dead_zone,
      cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 30, .01));

  // Repackage corner positions.
  std::vector<cv::Vec2f> points;
  points.reserve(_w * _h);
  for(size_t i = 0; i < _w * _h; ++i)
    points.emplace_back(corners.at<float>(i, 0), corners.at<float>(i, 1));

  return points;
}


std::vector<cv::Vec3f>
create_world_points(
    const size_t _w,
    const size_t _h,
    const double _square_size
) {
  std::vector<cv::Vec3f> target_points;
  target_points.reserve(_w * _h);

  for(size_t i = 0; i < _h; ++i)
    for(size_t j = 0; j < _w; ++j)
      target_points.emplace_back(j * _square_size, i * _square_size, 0.);

  return target_points;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~ camera_calibration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*------------------------------ Construction --------------------------------*/

camera_calibration::
camera_calibration()
{
  intrinsics = cv::Mat::eye(3, 3, CV_32F);
  distortion = cv::Mat::zeros(5, 1, CV_32F);
}


camera_calibration::
camera_calibration(
    const cv::Mat& _intrinsics,
    const cv::Mat& _distortion,
    const cv::Size& _image_size
)
  : intrinsics(_intrinsics), distortion(_distortion), image_size(_image_size)
{ }


camera_calibration::
camera_calibration(
    const std::string& _filename
) {
  read_xml(_filename);
}


camera_calibration::
camera_calibration(
    const std::vector<cv::Mat>& _images,
    const size_t _x,
    const size_t _y,
    const double _square_size
) {
  std::cout << "\nProcessing calibration images "
            << _images.front().size()
            << "...\n";

  // To calibrate, we need one point in both image and world coords for each
  // intersection of four squares.
  const size_t num_x = _x - 1,
               num_y = _y - 1;

  // Detect the square intersections in each image.
  std::vector<std::vector<cv::Vec2f>> image_points;
  for(size_t i = 0; i < _images.size(); ++i) {
    const auto& image = _images[i];
    try {
      auto points = extract_image_points(image, num_x, num_y);
      image_points.push_back(std::move(points));
    }
    catch(std::exception& _e) {
      std::cout << "Error in processing image " << i << ":\n\t" << _e.what()
                << "\nSkipping this image.\n"
                << std::endl;
    }
  }

  // Assert that we got at least two usable images.
  if(image_points.size() < 2) {
    std::cout << "\nAt least two usable images are required to calibrate."
              << std::endl;
    std::exit(0);
  }

  // Generate world points for the target.
  auto target_points = create_world_points(num_x, num_y, _square_size);

  // We need one set of these for each calibration image.
  std::vector<std::vector<cv::Vec3f>> world_points;
  for(size_t i = 0; i < image_points.size(); ++i)
    world_points.push_back(target_points);

  // Get the image size.
  image_size = _images.front().size();

  // Initialize calibration parameters.
  intrinsics = cv::Mat::eye(3, 3, CV_32F);
  distortion = cv::Mat::zeros(5, 1, CV_32F);

  // Make space for the extrinsics.
  std::vector<cv::Mat> rotations, translations;

  // Run OpenCV's image calibration routine.
  const double error = cv::calibrateCamera(world_points, image_points,
      image_size, intrinsics, distortion, rotations, translations);

  std::cout << "\nCamera intrinsics:\n" << intrinsics
            << "\nError: " << error << "\n"
            << std::endl;
}


void
camera_calibration::
compute_rectify_maps()
{
  remaps.clear();
  remaps.resize(4);

  cv::initUndistortRectifyMap(intrinsics, distortion, cv::Mat(), intrinsics,
      image_size, CV_32FC1, remaps[0], remaps[1]);
}

/*---------------------------------- I/O -------------------------------------*/

void
camera_calibration::
read_xml(
    const std::string& _filename
) {
  cv::FileStorage fs(_filename, cv::FileStorage::READ);
  fs["intrinsics"] >> intrinsics;
  fs["distortion"] >> distortion;

  // Read the size into a buffer to avoid problems with older OpenCV (v 2.4.5).
  std::vector<int> sizeBuffer;
  fs["image_size"] >> sizeBuffer;
  image_size = cv::Size(sizeBuffer[0], sizeBuffer[1]);

  fs.release();
}


void
camera_calibration::
write_xml(
    const std::string& _filename
) {
  cv::FileStorage fs(_filename, cv::FileStorage::WRITE);
  fs << "intrinsics" << intrinsics;
  fs << "distortion" << distortion;
  fs << "image_size" << image_size;
  fs.release();
}

/*----------------------------------------------------------------------------*/
