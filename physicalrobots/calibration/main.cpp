#include <string>

#include "calibration.h"
#include "image-capture.h"
#include "stereo-calibration.h"


int
main(const int _argc, const char* _argv[])
{
  const std::string instructions = "To create a good camera calibration, "
    "position a checkerboard target at various points in the field of view and "
    "take a calibration image with space. Try to get a uniform coverage of the "
    "view space from both a near (about 2ft) and a far (about 6ft) distance.";

#if 1
  // This version is for a single camera.

  // Ensure we got the right number of parameters.
  if(_argc != 5 and _argc != 6) {
    std::cerr << "usage: calibrate <camera index> <squares wide> <squares high> "
              << "<square size> [output xml filename]"
              << std::endl;
    std::exit(-1);
  }

  std::cout << instructions << std::endl;

  // Extract command-line parameters.
  const size_t camera_index = std::stoul(_argv[1]),
               squares_wide = std::stoul(_argv[2]),
               squares_high = std::stoul(_argv[3]);
  const double square_size  = std::stod(_argv[4]);

  // Capture a set of calibration images.
  std::vector<cv::Mat> calibration_images = capture_images(camera_index);

  // Perform auto calibration.
  camera_calibration cal(calibration_images, squares_wide, squares_high,
      square_size);

  // If output was requested, generate it now.
  const std::string output_file = _argc == 6 ? _argv[5] : "calibration.xml";
  cal.write_xml(output_file);

#else
  // This version is for a stereo rig.

  // Ensure we got the right number of parameters.
  if(_argc != 6 and _argc != 7) {
    std::cerr << "usage: calibrate <camera 1 index> <camera 2 index> "
              << "<squares wide> <squares high> "
              << "<square size> [output xml filename]"
              << std::endl;
    std::exit(-1);
  }

  std::cout << instructions << std::endl;

  // Extract command-line parameters.
  const size_t camera_1_index = std::stoul(_argv[1]),
               camera_2_index = std::stoul(_argv[2]),
               squares_wide = std::stoul(_argv[3]),
               squares_high = std::stoul(_argv[4]);
  const double square_size  = std::stod(_argv[5]);

  // Calibrate each camera.

  auto images_1 = capture_images(camera_1_index);
  auto images_2 = capture_images(camera_2_index);

  auto cal1 = camera_calibration(images_1, squares_wide, squares_high,
      square_size);
  auto cal2 = camera_calibration(images_2, squares_wide, squares_high,
      square_size);

  cal1.write_xml("cam1.xml");
  cal2.write_xml("cam2.xml");

  // Calibrate the rig.

  // Capture a set of calibration images.
  auto calibration_images = capture_images({camera_1_index, camera_2_index});

  // Perform auto calibration.
  stereo_calibration cal(cal1, cal2,
      calibration_images[0], calibration_images[1],
      squares_wide, squares_high, square_size);

  cal.compute_rectify_maps();
  cal.write_xml("rig.xml");
#endif
}
