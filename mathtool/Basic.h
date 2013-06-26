#ifndef BASIC_H_
#define BASIC_H_

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

namespace mathtool {

  //constants
#ifndef PI
#define PI 3.1415926535897
#endif

#ifndef TWOPI
#define TWOPI 6.2831853071794
#endif

  // Return the square of a.
  template<typename T>
    inline T sqr(const T& a) {return a*a;}

  // Return sign of x 
  inline int sign(double x) {return x >=0 ? 1: -1;}

  // Angle conversions
  inline double degToRad(double x) {return x*180/PI;}

  inline double radToDeg(double x) {return x*PI/180;}

  // computes sqrt(a^2 + b^2) without destructive underflow or overflow
  // prerequisite: a and b are positive real numbers
  inline double pythag(double a, double b) {return a * std::sqrt(1. + sqr(b/a));}
}

#endif
