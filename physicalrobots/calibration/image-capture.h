#ifndef _IMAGE_CAPTURE_H_
#define _IMAGE_CAPTURE_H_

#include <cstddef>
#include <vector>

#include <opencv2/opencv.hpp>


/// Capture a series of images from a designated camera.
/// @param _camera_index The index of the camera to use.
/// @param _maps Remaps for the camera (or empty for none).
/// @return A vector of the captured images.
std::vector<cv::Mat>
capture_images(
    const size_t _camera_index,
    const std::vector<cv::Mat>& _map = {}
);


/// Capture a series of images simultaneously from several cameras.
/// @param _camera_indexs The indexes of the cameras to use.
/// @param _maps Remaps for each camera (or empty for none).
/// @return A vector of the captured images for each camera.
std::vector<std::vector<cv::Mat>>
capture_images(
    const std::vector<size_t>& _camera_indexes,
    const std::vector<cv::Mat>& _maps = {}
);

#endif
