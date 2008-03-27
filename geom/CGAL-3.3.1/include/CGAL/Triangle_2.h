// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.3-branch/Kernel_23/include/CGAL/Triangle_2.h $
// $Id: Triangle_2.h 35642 2006-12-27 23:26:06Z spion $
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_TRIANGLE_2_H
#define CGAL_TRIANGLE_2_H

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Triangle_2 : public R_::Kernel_base::Triangle_2
{
  typedef typename R_::Point_2          Point_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;
  typedef typename R_::Kernel_base::Triangle_2  RTriangle_2;

  typedef Triangle_2                            Self;
  BOOST_STATIC_ASSERT((boost::is_same<Self, typename R_::Triangle_2>::value));

public:

  typedef RTriangle_2 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_                          R;
  typedef typename R::FT               FT;

  Triangle_2() {}

  Triangle_2(const RTriangle_2& t)
      : RTriangle_2(t) {}

  Triangle_2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
      : RTriangle_2(typename R::Construct_triangle_2()(Return_base_tag(), p,q,r)) {}

  FT
  area() const
  {
    return R().compute_area_2_object()(vertex(0), vertex(1), vertex(2));
  }

  Orientation
  orientation() const
  {
    return R().orientation_2_object()(vertex(0), vertex(1), vertex(2));
  }

  Bounded_side
  bounded_side(const Point_2 &p) const
  {
    return R().bounded_side_2_object()(*this,p);
  }

  Oriented_side
  oriented_side(const Point_2 &p) const
  {
    return R().oriented_side_2_object()(*this,p);
  }

  bool
  operator==(const Triangle_2 &t) const
  {
    return R().equal_2_object()(*this,t);
  }

  bool
  operator!=(const Triangle_2 &t) const
  {
    return !(*this == t);
  }

  typename Qualified_result_of<typename R::Construct_vertex_2, Triangle_2, int>::type
  vertex(int i) const
  {
    return R().construct_vertex_2_object()(*this,i);
  }

  typename Qualified_result_of<typename R::Construct_vertex_2, Triangle_2, int>::type
  operator[](int i) const
  {
    return vertex(i);
  }

  bool
  has_on_bounded_side(const Point_2 &p) const
  {
    return bounded_side(p) == ON_BOUNDED_SIDE;
  }

  bool
  has_on_unbounded_side(const Point_2 &p) const
  {
    return bounded_side(p) == ON_UNBOUNDED_SIDE;
  }

  bool
  has_on_boundary(const Point_2 &p) const
  {
    return bounded_side(p) == ON_BOUNDARY;
  }

  bool
  has_on_negative_side(const Point_2 &p) const
  {
    return oriented_side(p) == ON_NEGATIVE_SIDE;
  }

  bool
  has_on_positive_side(const Point_2 &p) const
  {
    return oriented_side(p) == ON_POSITIVE_SIDE;
  }

  bool
  is_degenerate() const
  {
    return R().collinear_2_object()(vertex(0), vertex(1), vertex(2));
  }

  Bbox_2
  bbox() const
  {
    return R().construct_bbox_2_object()(*this);
  }

  Triangle_2
  opposite() const
  {
    return R().construct_opposite_triangle_2_object()(*this);
  }

  Triangle_2
  transform(const Aff_transformation_2 &t) const
  {
    return Triangle_2(t.transform(vertex(0)),
                      t.transform(vertex(1)),
                      t.transform(vertex(2)));
  }

};


template < class R >
std::ostream &
operator<<(std::ostream &os, const Triangle_2<R> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0] << t[1]  << t[2];
    default:
        return os<< "Triangle_2(" << t[0] << ", "
                 << t[1] << ", " << t[2] <<")";
    }
}

template < class R >
std::istream &
operator>>(std::istream &is, Triangle_2<R> &t)
{
    typename R::Point_2 p, q, r;

    is >> p >> q >> r;

    if (is)
        t = Triangle_2<R>(p, q, r);
    return is;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLE_2_H
