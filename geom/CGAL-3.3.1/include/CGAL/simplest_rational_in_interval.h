// Copyright (c) 2002  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.3-branch/Number_types/include/CGAL/simplest_rational_in_interval.h $
// $Id: simplest_rational_in_interval.h 37955 2007-04-05 13:02:19Z spion $
//
//
// Author(s)     : Andreas Fabri, Susan Hert, Sylvain Pion

#ifndef CGAL_SIMPLEST_RATIONAL_IN_INTERVAL_H
#define CGAL_SIMPLEST_RATIONAL_IN_INTERVAL_H


#include <CGAL/number_type_basic.h>
#include <CGAL/to_rational.h>
#include <cassert>
#include <climits>
#include <cmath>


CGAL_BEGIN_NAMESPACE

/* simplest_rational_in_interval(x,y) returns the rational number with
     the smallest denominator in the interval [x,y].  See Knuth,
     "Seminumerical algorithms", page 654, answer to exercise
     4.53-39. */

template <class Rational>
Rational
simplest_rational_in_interval(double x, double y) {

    typedef Fraction_traits<Rational> FT;
    typedef typename FT::Is_fraction Is_fraction;
    typedef typename FT::Numerator_type Numerator_type;
    typedef typename FT::Denominator_type Denominator_type;
    typedef typename FT::Decompose Decompose;
    typedef typename FT::Compose Compose;

    // Must be a fraction
    BOOST_STATIC_ASSERT((::boost::is_same<Is_fraction, Tag_true>::value));
    // Numerator_type,Denominator_type must be the same
    BOOST_STATIC_ASSERT((::boost::is_same<Numerator_type, Denominator_type>::value));


  if(x == y){
    return to_rational<Rational>(x);
  }

  if(x > y){
    std::swap(x,y);
  }

  Rational r;  // Return value.
  Numerator_type r_numerator, r_denominator;
  // Deal with negative arguments.  We only have to deal with the case
  // where both x and y are negative -- when exactly one is negative
  // the best rational in the interval [x,y] is 0.
  if (x < 0 && y < 0) {
    // Both arguments are negative: solve positive case and negate
    return  - simplest_rational_in_interval<Rational>(std::fabs(x),std::fabs(y));
  } else if (x <= 0 || y <= 0) {
    // One argument is 0, or arguments are on opposite sides of 0:
    // simplest rational in interval is 0 exactly.
    r_numerator = 0;
    r_denominator = 1;
  } else { // x > 0 && y > 0
    double xc = CGAL_CLIB_STD::floor(1/x); // First coefficient of cf for x.
    double xr = CGAL_CLIB_STD::fmod(1/x,1); // Remaining fractional part of x.
    double yc = CGAL_CLIB_STD::floor(1/y); // First coefficient of cf for y.
    double yr = CGAL_CLIB_STD::fmod(1/y,1); // Remaining fractional part of y.

    if (xc < yc) {
      // Return 1/(xc+1).
      r_numerator = 1;
      r_denominator = xc + 1;
    } else if (yc < xc) {
      // Return 1/(yc+1).
      r_numerator = 1;
      r_denominator = yc + 1;
    } else  {  // xc == yc
      // Recurse to find s, the rational with the lowest denominator
      //  between xr and yr.
      Rational  s(simplest_rational_in_interval<Rational>(xr,yr));

      // Return 1/(xc + s).

      Numerator_type xc_rt(xc);
      Numerator_type s_num,s_den;
      Decompose()(s,s_num,s_den);
      r_numerator = s_den;
      r_denominator = s_num + xc_rt * s_den;
    }
  }

  return Compose()(r_numerator, r_denominator);
}

CGAL_END_NAMESPACE

#endif //CGAL_SIMPLEST_RATIONAL_IN_INTERVAL_H