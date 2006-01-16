#include <iostream.h>
#include <fstream.h>
#include <string.h>

#include "vclip.h"

int main(int argc, char **argv)
{
  Polyhedron *p1 = new Polyhedron;
  p1->addVertex("", Vect3( -1  ,  1  ,  1  ) );
  p1->addVertex("", Vect3( -1  ,  1  , -1  ) );
  p1->addVertex("", Vect3( -1  , -1  ,  1  ) );
  p1->addVertex("", Vect3( -1  , -1  , -1  ) );
  p1->addVertex("", Vect3(  1  ,  1  ,  1  ) );
  p1->addVertex("", Vect3(  1  ,  1  , -1  ) );
  p1->addVertex("", Vect3(  1  , -1  ,  1  ) );
  p1->addVertex("", Vect3(  1  , -1  , -1  ) );
  p1->buildHull();
  PolyTree *rob=new PolyTree;			// cube
  rob->setPoly(p1);				// ROBOT

  Polyhedron *p2 = new Polyhedron;
  p2->addVertex("", Vect3( -1  ,  1  ,  1  ) );
  p2->addVertex("", Vect3( -1  ,  1  , -1  ) );
  p2->addVertex("", Vect3( -1  , -1  ,  1  ) );
  p2->addVertex("", Vect3( -1  , -1  , -1  ) );
  p2->addVertex("", Vect3(  1  ,  1  ,  1  ) );
  p2->addVertex("", Vect3(  1  ,  1  , -1  ) );
  p2->addVertex("", Vect3(  1  , -1  ,  1  ) );
  p2->addVertex("", Vect3(  1  , -1  , -1  ) );
  p2->buildHull();
  PolyTree *obst=new PolyTree;			// cube
  obst->setPoly(p2);				// OBSTACLE


//----------------------------
// want to get from this to whatever vclip understands
//----------------------------

//  Cfg c;      //robot  <--- new!    pos&ori of bodyfixedframe robot
//  Cfg c_obst; //obst   <--- fixed   pos&ori of bodyfixedframe obst
//
//  Cfg Q  = c - c_obst; <--- feed this to vclip



  double Q[6]={0,0,5,   0,0,0.5};  // (xyz rpy)


//----------------------------
// this is what i wanted to do, or something simple
//----------------------------
  Vect3 XYZ(Q[0],Q[1],Q[2]);

  double const TWOPI = 3.14159*2.0;

  Quat roll (Q[3]*TWOPI,Vect3::I);
  Quat pitch(Q[4]*TWOPI,Vect3::J);
  Quat yaw  (Q[5]*TWOPI,Vect3::K);

  Quat combo;
  combo.mult(roll,pitch);
  combo.postmult(yaw);
  VclipPose Xbeauty(combo,XYZ);
  
//----------------------------
// this is what daniel does (although, I've probably made errors)
//----------------------------

  double c=Q[3]*TWOPI,
         b=Q[4]*TWOPI,
         a=Q[5]*TWOPI,

         ca=cos(a),
         cb=cos(b),
         cc=cos(c),
         sa=sin(a),
         sb=sin(b),
         sc=sin(c),

         tmp1=ca*sb,
         tmp2=sa*sb;

  Mat3 M;
       M.setXcol(Vect3( ca*cb, 
                        sa*cb, 
                          -sb) );
       M.setYcol(Vect3( tmp1*sc - sa*cc,
                        tmp2*sc + ca*cc,
                                  cb*sc) );
       M.setZcol(Vect3( tmp1*cc + sa*sc,
                        tmp2*cc - ca*sc,
                                  cb*cc) );

  VclipPose Xtrans(M,XYZ);

  Quat comparethis(M);

  cout<<"\n Matrix's \n";
  comparethis.print(cout);
  cout<<"\n\n Quaternion's \n";
  combo.print(cout);

  if (comparethis == combo) cout<<"\n\n methods are the same!\n";


  Vect3 cp1, cp2;   // closest points between bodies, in local frame
  ClosestFeaturesHT closestFeaturesHT(3000);
  Real dist = PolyTree::vclip(rob,obst,Xbeauty,closestFeaturesHT, cp1, cp2);

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

