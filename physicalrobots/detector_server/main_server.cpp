#define PI 3.1415926

#include <unistd.h>

#include "ArucoDetector.h"
#include <iostream>
#include <string>
#include <unistd.h>
#include <thread>

#include "packet.h"
#include "tcp_socket.h"

using namespace std;
using namespace nonstd;


// Run the detector locally (do not accept tcp connections, just show the
// markers and measured data). This is useful for test/debug of the detector
// server.
void
local(const std::string& _calibrationFile)
{
  // Create an aruco detector object.
  ArucoDetector detector(_calibrationFile);
  detector.Display(true);

  while(true) {
    // Detect markers and report the number seen.
    vector<MarkerPercept> percepts = detector();
    cout << "Spoted " << percepts.size() << " markers."
         << std::endl;

    // If we saw markers, send them to the client.
    for(const auto& p : percepts) {
      cout << "\tMarker " << p.ID()
           << " at distance " << p.Distance()
           << ", heading " << p.Heading()
           << ", angle " << p.Angle()
           << endl;
    }

    usleep(100000);
  }
}


// Run the detector as a server (respond to tcp connections and send marker
// data).
void
server(const std::string& _calibrationFile)
{
  tcp_socket server;

  // Create an aruco detector object.
  ArucoDetector detector(_calibrationFile);
  detector.Display(true);

  // Define the server handling function.
  tcp_socket::handler_function f = [&](tcp_socket&& _client) {
    while(true) {
      // Wait for the client to signal detection.
      std::cout << "Waiting on client signal..."
                << std::endl;
      bool go = false;
      _client >> go;

      if(!go)
        break;

      // Detect markers and report the number seen.
      vector<MarkerPercept> percepts = detector();
      cout << "Spoted " << percepts.size() << " markers."
           << std::endl;
      _client << (char)(percepts.size());

      // If we saw markers, send them to the client.
      for(const auto& p : percepts) {
        cout << "\tMarker " << p
             << endl;

        auto local = p.GetCameraPositionInMarkerFrame();

        _client << packet(p.ID(), packet::s_factor * local[0],
                                  packet::s_factor * local[1],
                                  packet::s_factor * local[2]);
      }
    }
  };

  server.set_handler(f);
  server.listen("4002", 1000, false);
}


int main(int _argc, char* _argv[])
{
  if(_argc != 2)
  {
    std::cout << "usage: ./detector <calibration file>" << std::endl;
    exit(1);
  }

  const std::string calibrationFile = _argv[1];
  //local(calibrationFile);
  server(calibrationFile);
}
