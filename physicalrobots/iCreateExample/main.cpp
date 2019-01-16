#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <thread>

#include <unistd.h>

#include "RobotController.h"

using namespace std;


int main(int argc, char *argv[])
{
  if(argc != 2)
  {
     cerr << "usage: ./example <robot_ip>" << endl;
     exit(-1);
  }


  RobotController r(argv[1]);
  r.StartCommandQueue();
  r.EnqueueCommand(.5,  0, 2);
  r.EnqueueCommand( 0, .5, 2);
  r.EnqueueCommand(.5,  0, 2);
  r.EnqueueCommand(0, 0, .1);
  usleep(10 * 1e6);
}
