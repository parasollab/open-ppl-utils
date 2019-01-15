// -*- c++ -*-
//
// Copyright (c) 2000-2009, Texas Engineering Experiment Station (TEES), a
// component of the Texas A&M University System.
//
// All rights reserved.
//
// The information and source code contained herein is the exclusive
// property of TEES and may not be disclosed, examined or reproduced
// in whole or in part without explicit written authorization from TEES.
//
// Copyright 1997, 1998, 1999 University of Notre Dame.
// Authors: Andrew Lumsdaine, Jeremy G. Siek, Lie-Quan Lee
//
// This file is part of the Matrix Template Library
//
// You should have received a copy of the License Agreement for the
// Matrix Template Library along with the software;  see the
// file LICENSE.  If not, contact Office of Research, University of Notre
// Dame, Notre Dame, IN  46556.
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//
//===========================================================================

#include "mtl/utils.h"
#include "mtl/matrix.h"
#include "matrix_test.h"

// matrix_attr.h is generated by make_and_test.pl and defines
// NUMTYPE, SHAPE, STORAGE, ORIEN, and TESTNAME
// you can create your own for testing purposes
#include "matrix_attr.h"


template <class Matrix>
inline bool const_oned_test(const Matrix& A, std::string test_name)
{
  typedef typename mtl::matrix_traits<Matrix>::value_type T;
  typedef typename mtl::matrix_traits<Matrix>::size_type Int;
  {
    /* read it out */
    T c = T(0);
    for (Int i = 0; i < A.noneds(); ++i) {
      for (typename Matrix::OneDRef::const_iterator j = A[i].begin();
	   j != A[i].end(); ++j) {
	c = c + T(1);
	if (*j != c) {
	  std::cerr << "**** FAILED: (const oned) "
	       << test_name.c_str() << " ****" << std::endl;
	  return false;
	}
      }
    }
  }
  return true;
}

template <class Matrix>
inline bool oned_test(Matrix& A, std::string test_name)
{
  typedef typename mtl::matrix_traits<Matrix>::value_type T;
  typedef typename mtl::matrix_traits<Matrix>::size_type Int;

  {
    /* read it out */
    T c = T(0);
    for (Int i = 0; i < A.noneds(); ++i) {
      for (typename Matrix::OneDRef::iterator j = A[i].begin();
	   j != A[i].end(); ++j) {
	c = c + T(1);
	if (*j != c) {
	  std::cerr << "**** FAILED: (oned) "
	       << test_name.c_str() << " ****" << std::endl;
	  return false;
	}
      }
    }
  }
  bool ret = const_oned_test(A, test_name);
  if (ret) {
    std::cout << test_name.c_str() << " passed oned_test" << std::endl;
    return true;
  } else
    return false;
}


template <class Matrix>
void
do_test(Matrix& A, std::string test_name)
{
  using namespace mtl;

  typedef typename mtl::matrix_traits<Matrix>::value_type T;
  typedef typename mtl::matrix_traits<Matrix>::size_type Int;

  iterator_fill(A);

  oned_test(A, test_name);
}


int
main(int argc, char* argv[])
{
  if (argc < 5) {
    std::cerr << "matrix_test <M> <N> <SUB> <SUPER>" << std::endl;
    return -1;
  }

  using namespace mtl;
  using std::string;

  const int M = atoi(argv[1]);
  const int N = atoi(argv[2]);
  const int SUB = atoi(argv[3]);
  const int SUP = atoi(argv[4]);

  std::cout << "M: " << M << " N: " << N 
       << " SUB: " << SUB << " SUPER: " << SUP << std::endl;

  typedef matrix<NUMTYPE, SHAPE, STORAGE, ORIEN>::type Matrix;

  string test_name = TESTNAME;
  Matrix* a = 0;

  create_and_run(M, N, SUB, SUP, test_name, a, Matrix::shape());

  return 0;
}
