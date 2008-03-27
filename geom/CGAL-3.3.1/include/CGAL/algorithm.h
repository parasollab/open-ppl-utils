// Copyright (c) 2003  Utrecht University (The Netherlands),
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.3-branch/STL_Extension/include/CGAL/algorithm.h $
// $Id: algorithm.h 37883 2007-04-03 16:07:32Z ameyer $
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>

#ifndef CGAL_ALGORITHM_H
#define CGAL_ALGORITHM_H 1

#include <CGAL/basic.h>
#include <CGAL/copy_n.h>
#include <algorithm>

#include <iosfwd>

CGAL_BEGIN_NAMESPACE

template <class ForwardIterator>
inline
ForwardIterator
successor( ForwardIterator it )
{
  return ++it;
}

template <class BidirectionalIterator>
inline
BidirectionalIterator
predecessor( BidirectionalIterator it )
{
  return --it;
}

template < class ForwardIterator >
std::pair< ForwardIterator, ForwardIterator >
min_max_element(ForwardIterator first, ForwardIterator last)
{
  typedef std::pair< ForwardIterator, ForwardIterator > FP;
  FP result(first, first);
  if (first != last)
    while (++first != last) {
      if (*first < *result.first)
        result.first = first;
      if (*result.second < *first)
        result.second = first;
    }
  return result;
}

template < class ForwardIterator, class CompareMin, class CompareMax >
std::pair< ForwardIterator, ForwardIterator >
min_max_element(ForwardIterator first,
                ForwardIterator last,
                CompareMin comp_min,
                CompareMax comp_max)
{
  typedef std::pair< ForwardIterator, ForwardIterator > FP;
  FP result(first, first);
  if (first != last)
    while (++first != last) {
      if (comp_min(*first, *result.first))
        result.first = first;
      if (comp_max(*result.second, *first))
        result.second = first;
    }
  return result;
}

template < class ForwardIterator, class Predicate >
ForwardIterator
min_element_if(ForwardIterator first,
               ForwardIterator last,
               Predicate pred)
{
  ForwardIterator result = first = std::find_if(first, last, pred);
  if (first != last)
    while (++first != last)
      if (*first < *result && pred(*first))
        result = first;
  return result;
}

template < class ForwardIterator, class Compare, class Predicate >
ForwardIterator
min_element_if(ForwardIterator first,
               ForwardIterator last,
               Compare comp,
               Predicate pred)
{
  ForwardIterator result = first = std::find_if(first, last, pred);
  if (first != last)
    while (++first != last)
      if (comp(*first, *result) && pred(*first))
        result = first;
  return result;
}

template < class ForwardIterator, class Predicate >
ForwardIterator
max_element_if(ForwardIterator first,
               ForwardIterator last,
               Predicate pred)
{
  ForwardIterator result = first = std::find_if(first, last, pred);
  if (first != last)
    while (++first != last)
      if (*result < *first && pred(*first))
        result = first;
  return result;
}

template < class ForwardIterator, class Compare, class Predicate >
ForwardIterator
max_element_if(ForwardIterator first,
               ForwardIterator last,
               Compare comp,
               Predicate pred)
{
  ForwardIterator result = first = std::find_if(first, last, pred);
  if (first != last)
    while (++first != last)
      if (comp(*result, *first) && pred(*first))
        result = first;
  return result;
}



/*! \brief lexicographic comparison of the two ranges using the \a cmp
    function object.

    Compares the two ranges \c [first1,last1) and \c [first2,last2)
    lexicographically and returns one of the \c CGAL::Comparison_result enum
    values respectively:
      - \c CGAL::SMALLER
      - \c CGAL::EQUAL
      - \c CGAL::LARGER

    \pre The \c value_type of \a InputIterator1 must be convertible
    to the \c first_argument_type of \c BinaryFunction.
    The \c value_type of \a InputIterator2 must be convertible
    to the \c second_argument_type of \c BinaryFunction.
    The \c result_type of \c BinaryFunction must be convertible to
    \c CGAL::Comparison_result.
 */

template <class InputIterator1, class InputIterator2, class BinaryFunction>
CGAL::Comparison_result
lexicographical_compare_three_valued( InputIterator1 first1, InputIterator1 last1,
             InputIterator2 first2, InputIterator2 last2,
             BinaryFunction cmp) {
    while ( first1 != last1 && first2 != last2) {
        CGAL::Comparison_result result = cmp( *first1, *first2);
        if ( result != CGAL::EQUAL)
            return result;
        ++first1;
        ++first2;
    }
    if ( first1 != last1)
        return CGAL::LARGER;
    if ( first2 != last2)
        return CGAL::SMALLER;
    return CGAL::EQUAL;
}

/*! \brief output iterator range to a stream, with separators

    The iterator range \c [first,beyond) is written
    to \c os (obeying CGAL I/O modes). Each element is bracketed by
    \c pre and \c post (default: empty string). Adjacent values are
    spearated by \c sep (default: ", ", i.e. comma space).
    The stream \c os is returned in its new state after output.

    Example:
<PRE>
    int a[] = {1, 2, 3};
    output_range(std::cout, a, a+3, ":", "(", ")");
</PRE>
    produces \c (1):(2):(3)
 */
template <class InputIterator>
std::ostream& 
output_range(std::ostream& os,
             InputIterator first, InputIterator beyond,
             const char* sep = ", ", const char* pre = "", const char* post = "") 
{
    InputIterator it = first;
    if (it != beyond) {
        os << pre << oformat(*it) << post;
        while (++it != beyond) os << sep << pre << oformat(*it) << post;
    }
    return os;
}



CGAL_END_NAMESPACE

#endif // CGAL_ALGORITHM_H //
// EOF //
