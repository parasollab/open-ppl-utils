#include <iostream.h>
#include <fstream.h>
#include <string.h>

#include "vclip.h"

int main(int argc, char **argv)
{
  Polyhedron *p1 = new Polyhedron;
  p1->addVertex("", Vect3(0, 0, 0) );
  p1->addVertex("", Vect3(1, 0, 0) );
  p1->addVertex("", Vect3(0, 1, 0) );
  p1->addVertex("", Vect3(0, 0, 1) );
  p1->buildHull();


  PolyTree *rob=new PolyTree;
  rob->setPoly(p1);

  Polyhedron *p2 = new Polyhedron;
  p2->addVertex("", Vect3(-0.5, 0.5, 0.5) );
  p2->addVertex("", Vect3(-0.5, 0.5,-0.5) );
  p2->addVertex("", Vect3(-0.5,-0.5, 0.5) );
  p2->addVertex("", Vect3(-0.5,-0.5,-0.5) );
  p2->addVertex("", Vect3( 0.5, 0.5, 0.5) );
  p2->addVertex("", Vect3( 0.5, 0.5,-0.5) );
  p2->addVertex("", Vect3( 0.5,-0.5, 0.5) );
  p2->addVertex("", Vect3( 0.5,-0.5,-0.5) );
  p2->buildHull();
  PolyTree *obst=new PolyTree;
  obst->setPoly(p2);

  VclipPose X12= VclipPose::ID;    // xform from body1 frame to body2 frame
                                        // local & global frames coincide here..
  Vect3 cp1, cp2;   // closest points between bodies, in local frame
  ClosestFeaturesHT closestFeaturesHT(3000);
  Real dist = PolyTree::vclip(rob,obst,X12,closestFeaturesHT, cp1, cp2);

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

