// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.3-branch/Skin_surface_3/include/CGAL/Skin_surface_polyhedral_items_with_face_information.h $
// $Id: Skin_surface_polyhedral_items_with_face_information.h 28711 2006-02-23 10:59:28Z nicokruithof $
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SKIN_SURFACE_POLYHEDRAL_ITEMS_WITH_FACE_INFORMATION_H
#define CGAL_SKIN_SURFACE_POLYHEDRAL_ITEMS_WITH_FACE_INFORMATION_H

CGAL_BEGIN_NAMESPACE

template <class Refs, class TriangulatedMixedComplex3>
struct Skin_Surface_polyhedral_face : public CGAL::HalfedgeDS_face_base<Refs> {
  typedef typename TriangulatedMixedComplex3::Cell_handle Triang_Cell_handle;

  Triang_Cell_handle triang_ch;
};

template < class TriangulatedMixedComplex3 >
class Skin_surface_polyhedral_items_with_face_information_3
  : public Polyhedron_items_3 {
  
  template <class Refs, class Traits>
  struct Face_wrapper {
    typedef Skin_Surface_polyhedral_face<Refs, TriangulatedMixedComplex3> Face;
  };
};

CGAL_END_NAMESPACE

#endif // CGAL_SKIN_SURFACE_POLYHEDRAL_ITEMS_WITH_FACE_INFORMATION_H
