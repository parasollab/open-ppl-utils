#include "image-capture.h"

#include <iostream>

#include "calibration.h"


std::vector<cv::Mat>
capture_images(
    const size_t _camera_index,
    const std::vector<cv::Mat>& _maps
) {
  auto images = capture_images(std::vector<size_t>{_camera_index}, _maps);
  return images.front();
}


std::vector<std::vector<cv::Mat>>
capture_images(
    const std::vector<size_t>& _camera_indexes,
    const std::vector<cv::Mat>& _maps
)
{
  // We will capture images from this many cameras at the same time.
  const size_t num = _camera_indexes.size();

  // Create video capture object for each camera.
  std::vector<cv::VideoCapture> capture(num);
  for(size_t i = 0; i < num; ++i) {
    capture[i].open(_camera_indexes[i]);

    if(!capture[i].isOpened()) {
      std::cerr << "Could not open video stream for camera " << _camera_indexes[i]
                << "." << std::endl;
      std::exit(-1);
    }
  }

  // For each camera, create a vector of images to store the calibration shots.
  std::vector<cv::Mat> frames(num);
  std::vector<std::vector<cv::Mat>> images(num);

  // Capture as many images as the user will provide.
  std::cout << "Press *Space* to capture an image, or *q* to end input.\n";

  bool go = true;
  const bool remaps = _maps.size();
  while(go) {
    // Capture the current frames.
    cv::Mat buffer;
    for(size_t i = 0; i < num; ++i) {
      capture[i] >> buffer;

      // Remap if we are using maps.
      if(remaps) {
        cv::remap(buffer, frames[i], _maps[2*i], _maps[2*i + 1],
            cv::INTER_LINEAR);
        const cv::Size size = frames[i].size();
        for(int j = size.height / 10 - 1; j < size.height;
            j += size.height / 10)
          cv::line(frames[i], cv::Point(0, j), cv::Point(size.width, j),
              cv::Scalar(0, 255, 0));
      }
      else
        frames[i] = buffer.clone();

      cv::imshow("Camera " + std::to_string(_camera_indexes[i]) + " Feed",
          frames[i]);
    }

    // Delay for 30 ms before showing the next frame. Exit on key press.
    const char key = cv::waitKey(30);
    switch(key) {
      case 'q':
      case 'Q':
        go = false;
        std::cout << "\tInput ended.\n";
        break;
      case ' ':
        for(size_t i = 0; i < num; ++i) {
          // Capture picture.
          images[i].push_back(frames[i].clone());

          // Try to detect and draw points.
          try {
            auto pts = extract_image_points(frames[i], 6, 9);
            for(const auto& p : pts)
              cv::circle(frames[i], cv::Point(int(p[0]), int(p[1])), 4,
                  cv::Scalar(0,0,255));
          }
          catch(std::exception& _e) {
            std::cerr << "\n" << _e.what() << "\n";
          }

          // Show modified picture.
          cv::imshow("Camera " + std::to_string(_camera_indexes[i]) + " Feed",
              frames[i]);
        }
        std::cout << "\t" << images.front().size() << " images captured.\n";
        cv::waitKey(0);
        break;
      default:;
    }
  }

  std::cout << "Captured " << images.front().size() << " images.\n";

  return images;
}
