#include <iostream.h>
#include <fstream.h>
#include <string.h>

#include "vclip.h"

int main(int argc, char **argv)
{
  PolyTreeLibrary polyTreeLibrary;
  loadPolyTreeFile("simple.pt", polyTreeLibrary);

  ClosestFeaturesHT closestFeaturesHT(3000); 

  PolyTree *b1=polyTreeLibrary.create("tetra"),
           *b2=polyTreeLibrary.create("unit-cube");

  VclipPose X12= VclipPose::ID;    // xform from body1 frame to body2 frame
					// local & global frames coincide here...

  Real dist;        // distance between bodies
  Vect3 cp1, cp2;   // closest points between bodies, in local frame

  dist = PolyTree::vclip(b1,b2,X12,closestFeaturesHT, cp1, cp2);

  cout<< "\n\tdist = " << dist;
  if (dist > 0.0)
	cout<< "\n\tNot in collision\n\n";
	// this dist is a lowerbound to real distance between these bodies
	// witness pairs may not be on real surfaces but rather on 
	// chull approx's
  else
	cout<< "\n\tCollision!!!\n\n";

  // in theory, equal 0.0 => contact

}

