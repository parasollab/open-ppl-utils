What is in this repo:
- aruco: The aruco detector library, used for marker detection.
- player: The player library, used to control the iCreate robots.
- detector_server: The aruco detector server which should run on a netbook (uses
  aruco). It will respond to TCP requests by detecting markers and sending their
  relative distance and heading back to the client.
- iCreateExample: An example controller for the iCreate using the player
  library. PMPL already has an iCreate controller, but this one is separate and
  easier to use for testing the robot itself.
- create.cfg: A player config file for the iCreate.

To run the iCreate:
- Start the player client on the netbook with 'sudo player create.cfg'. This
  creates a tcp interface for commanding the robot through a client controller
  like the one in pmpl or iCreateExample.
- Separately start the aruco detector server if you want to detect markers. This
  will be the executable called detector_server/detector; it also creates a TCP
  interface for getting marker data, which can be used by PMPL's
  ArucoDetectorInterface.

Note about building with old packages:
- If building the 'calibration' and/or 'detector_server' tools on a system with
  ancient packages like CentOS, you will need to adjust the Makefiles by
  removing the -lopencv_videoio entry in the OPENCV_LIBS variable (this library
  does not exit in opencv2).
