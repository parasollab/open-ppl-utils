// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.3-branch/Minkowski_sum_2/include/CGAL/Minkowski_sum_2/Offset_conv_2.h $
// $Id: Offset_conv_2.h 37897 2007-04-03 18:34:02Z efif $
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_OFFSET_CONV_H
#define CGAL_OFFSET_CONV_H

#include <CGAL/Minkowski_sum_2/Union_of_curve_cycles_2.h>
#include <list>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class for computing the offset of a given polygon by a given radius
 * by constructing a single convolution cycle and computing its interior.
 */
template <class Base_>
class Offset_by_convolution_2 : private Base_
{
private:

  typedef Base_                                          Base;

public:

  typedef typename Base::Basic_kernel                    Kernel;
  typedef typename Base::Basic_NT                        NT;
  typedef typename Base::Polygon_2                       Polygon_2;  
  typedef typename Base::Offset_polygon_2                Offset_polygon_2;
  
private:

  typedef typename Base::Labeled_traits_2                Labeled_traits_2; 
  typedef typename Base::Labeled_curve_2                 Labeled_curve_2;
  typedef std::list<Labeled_curve_2>                     Curves_list;

  typedef Union_of_curve_cycles_2<Labeled_traits_2,
                                  Offset_polygon_2>      Union_2;

public:

  /*! Constructor. */
  Offset_by_convolution_2 (const Base_& base) :
    Base (base)
  {}    

  /*!
   * Compute the offset of a simple polygon by a given radius.
   * Note that as the input polygon may not be convex, its offset may not be 
   * simply connected. The result is therefore represented as the outer
   * boundary of the Minkowski sum (which is always a simple offset polygon)
   * and a container of offset polygons, representing the holes in this "outer"
   * polygon.
   * \param traits Arrangement traits that can deal with line segments and 
   *               circular arcs.
   * \param pgn The polygon.
   * \param r The offset radius.
   * \param off_bound Output: The outer boundary of the offset polygon.
   * \param off_holes Output: An output iterator for the holes in the offset.
   * \pre The polygon is simple.
   * \return A past-the-end iterator for the holes container.
   */
  template <class OutputIterator>
  OutputIterator operator() (const Polygon_2& pgn,
                             const NT& r,
                             Offset_polygon_2& off_bound,
                             OutputIterator off_holes) const
  {
    CGAL_precondition (pgn.is_simple());

    // Compute the curves that form the single convolution cycle for the
    // given polygon.
    Curves_list                     cycle;

    _offset_polygon (pgn, r,
                     1,                       // The ID of the single cycle.
                     std::back_inserter (cycle));

    // Compute the union of the cycles that represent the offset polygon.
    Union_2     unite;

    off_holes = unite (cycle.begin(), cycle.end(),
                       off_bound, off_holes);

    return (off_holes);
  }

};


CGAL_END_NAMESPACE

#endif
