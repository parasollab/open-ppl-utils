// Copyright (c) 1997-2002  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.3-branch/Triangulation_2/demo/Triangulation_2/triangulation_2_edit_vertex.cpp $
// $Id: triangulation_2_edit_vertex.cpp 37003 2007-03-10 16:55:12Z spion $
//
//
// Author(s)     : Laurent Rineau

#include <CGAL/basic.h>

#ifdef CGAL_USE_QT

#include "triangulation_2_edit_vertex.h"

void triangulation_2_edit_vertex_helper::delete_vertex()
{
  delete_vertexi();
  emit(triangulation_changed());
};
void triangulation_2_edit_vertex_helper::move_vertex() { move_vertexi(); };
void triangulation_2_edit_vertex_helper::change_weight() { change_weighti(); };
void triangulation_2_edit_vertex_helper::stateChanged(int i){
  if(i==2)
    activate();
  else if(i == 0)
    deactivate();
}

#include "triangulation_2_edit_vertex.moc"

#endif