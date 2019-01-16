#define PI 3.1415926

#include <iostream>
#include <string>

#include "packet.h"
#include "tcp_socket.h"

#include <unistd.h>

using namespace std;
using namespace nonstd;


int
main(int _argc, char* _argv[])
{
  // Require exactly two arguments.
  assert_msg(_argc == 2, "usage: ./client <server_ip>");

  // Connect to the server.
  const std::string ip   = _argv[1],
                    port = "4002";
  tcp_socket server(ip, port);

  // Loop until Ctrl-C.
  while(true) {
    // Signal the server to detect.
    bool go = true;
    server << go;

    // Get the number of markers seen.
    char count = 0;
    server >> count;
    std::cout << "The detector saw " << (size_t)count << " markers.\n";

    // Get each marker.
    for(size_t i = 0; i < (size_t)count; i++) {
      packet p;
      server >> p;

      std::cout << "\nCamera position in marker frame:"
                << "\nX: " << p.x / packet::s_factor
                << "\nY: " << p.y / packet::s_factor
                << "\nT: " << p.t / packet::s_factor
                << std::endl;
    }

    usleep(100000);
  }
}
