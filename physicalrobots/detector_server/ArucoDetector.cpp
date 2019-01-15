#include "ArucoDetector.h"

#include "calibration.h"

#include <iostream>
#include <fstream>

using namespace std;


/*----------------------------- Construction ---------------------------------*/

ArucoDetector::
ArucoDetector(const string& _filename) {
  // Open and validate capture stream.
  m_capture.set(CV_CAP_PROP_FRAME_COUNT, 1);
  m_capture.open(0);
  if(!m_capture.isOpened()) {
    cerr << "ArucoDetector error: could not open video stream." << endl;
    exit(-1);
  }

  // Read calibration file if provided.
  if(_filename != "")
    ReadCalibrationFile(_filename);

  // Read an image from the webcam to set the camera image size.
  CaptureVideo();
  if(m_cameraParams.isValid())
    m_cameraParams.resize(m_image.size());

  // Set the marker dictionary to the default (ARUCO dictionary).
  m_detector.setDictionary(aruco::Dictionary::getTypeFromString("ARUCO"));

  // Set the detector thresholds. TODO: figure out what this means.
  m_detector.setThresholdParams(9., 10.);

  cout << "ArucoDetector initialized." << endl
       << "\tImage size: " << m_image.size() << endl;
}

/*---------------------------- Interface -------------------------------------*/

vector<MarkerPercept>
ArucoDetector::
operator()() {
  vector<MarkerPercept> output;

  auto markers = GetVisibleMarkers();
  for(const auto& m : markers)
    output.emplace_back(m);

  if(m_display)
    UpdateDisplay(markers);

  return output;
}


void
ArucoDetector::
Display(const bool _d) {
  m_display = _d;
}

/*---------------------------- Capture Helpers -------------------------------*/

void
ArucoDetector::
CaptureVideo(const size_t _skipFrames) {
  /// @TODO Remove _skipFrames and grab calls when we move to a newer OpenCV
  ///       version which supports setting the buffer size with
  ///       m_capture.set(CV_CAP_PROP_BUFFERSIZE, 1).
  try {
    // Grab enough frames to clear the buffer.
    for(size_t i = 0; i < _skipFrames; ++i)
      m_capture.grab();

    // Retrieve the last frame.
    m_capture >> m_image;
  }
  catch(cv::Exception& _e) {
    cerr << "ArucoDetector::CaptureVideo caught OpenCV exception: "
         << _e.what() << endl;
  }
}


vector<aruco::Marker>
ArucoDetector::
GetVisibleMarkers() {
  ///@todo This function currently only accomodates two fixed marker sizes. We
  ///      need to generalize it to handle multiple sizes.
  CaptureVideo();

  // The small marker size, like those on the robots (3.25").
  static constexpr float small = .0825;
  // The large marker size, like those on the walls of 407.
  static constexpr float big = .0825*2;

  vector<aruco::Marker> bigMarkers = m_detector.detect(m_image, m_cameraParams,
      big);

  /// Marker IDs 300 and 400 refer to small markers.
  for(auto iter = bigMarkers.begin(); iter != bigMarkers.end();){
    if(iter->id > 300 and iter->id < 400)
      iter = bigMarkers.erase(iter);
    else
      iter++;
  }

  vector<aruco::Marker> smallMarkers = m_detector.detect(m_image, m_cameraParams,
      small);
  for(auto iter = smallMarkers.begin(); iter != smallMarkers.end();){
    if(iter->id > 300 and iter->id < 400)
      iter++;
    else
      iter = smallMarkers.erase(iter);
  }

  bigMarkers.insert(bigMarkers.end(), smallMarkers.begin(), smallMarkers.end());
  return bigMarkers;
}

/*---------------------------- Display Helpers -------------------------------*/

void
ArucoDetector::
UpdateDisplay(const vector<aruco::Marker>& _markers) {
  try {
    DrawMarkers(_markers, m_image);
    cv::imshow("input", m_image);
    cv::imshow("thresholded", m_detector.getThresholdedImage());
    cv::waitKey(10);
  }
  catch(cv::Exception& _e) {
    cerr << "ArucoDetector::ShowDisplay caught OpenCV exception: "
         << _e.what() << endl;
  }
}


void
ArucoDetector::
DrawMarkers(const vector<aruco::Marker>& _markers, cv::Mat& _image) {
  // Draw each marker ID and outline.
  for(auto& m : _markers)
    m.draw(_image, cv::Scalar(0, 0, 255), 1);

  // Draw a 3D box and axis on each marker.
  if(m_cameraParams.isValid()){
    for(auto& m: _markers) {
      aruco::CvDrawingUtils::draw3dCube(_image, const_cast<aruco::Marker&>(m),
          m_cameraParams);
      aruco::CvDrawingUtils::draw3dAxis(_image, const_cast<aruco::Marker&>(m),
          m_cameraParams);
    }
  }
}

/*------------------------------ Setup Helpers -------------------------------*/

void
ArucoDetector::
ReadCalibrationFile(const std::string& _filename) {
  camera_calibration cal(_filename);
  m_cameraParams.setParams(cal.intrinsics, cal.distortion, cal.image_size);
}

/*----------------------------------------------------------------------------*/
