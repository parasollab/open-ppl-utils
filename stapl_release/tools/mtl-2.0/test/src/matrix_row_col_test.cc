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
#include "iter_ij_test.h"

// matrix_attr.h is generated by make_and_test.pl and defines
// NUMTYPE, SHAPE, STORAGE, ORIEN, and TESTNAME
// you can create your own for testing purposes
#include "matrix_attr.h"


template <class Matrix, class Shape>
bool rows_test(const Matrix& A, std::string test_name, Shape, row_tag)
{
  return iterator_operator_ij_test(A, test_name);
}

template <class Matrix>
bool rows_test(const Matrix& A, std::string test_name,
	       rectangle_tag, column_tag)
{
  return strided_iterator_operator_ij_test(A, test_name);
}

template <class Matrix, class Shape>
bool columns_test(const Matrix& A, std::string test_name, Shape, row_tag)
{
  return strided_iterator_operator_ij_test(A, test_name);
}

template <class Matrix>
bool columns_test(const Matrix& A, std::string test_name,
		  rectangle_tag, column_tag)
{
  return iterator_operator_ij_test(A, test_name);
}


template <class Matrix>
bool rows_columns_test(Matrix& A, std::string test_name, strideable)
{
  using mtl::rows;
  using mtl::columns;
  typedef typename mtl::matrix_traits<Matrix>::shape Shape;
  typedef typename mtl::matrix_traits<Matrix>::orientation Orien;
  bool ret1 = mtl::matrix_equal(rows(A), columns(A));
  bool ret2 = rows_test(rows(A), test_name, Shape(), Orien());
  bool ret3 = columns_test(columns(A), test_name, Shape(), Orien());
  if (ret1 && ret2 && ret3)
    std::cout << test_name.c_str() << " passed rows & columns test" << std::endl;
  return ret1 && ret2 && ret3;
}

template <class Matrix>
bool rows_columns_test(Matrix&, std::string test_name, not_strideable)
{
  std::cout << test_name.c_str() << " skipping rows & columns test" << std::endl;
  return true;
}

template <class Matrix>
bool rows_columns_test(Matrix& A, std::string test_name)
{
  typedef typename mtl::matrix_traits<Matrix>::strideability Stridable;
  return rows_columns_test(A, test_name, Stridable());
}


template <class Matrix>
void
do_test(Matrix& A, std::string test_name)
{
  using namespace mtl;

  typedef typename mtl::matrix_traits<Matrix>::value_type T;
  typedef typename mtl::matrix_traits<Matrix>::size_type Int;

  matrix_fill(A);

  rows_columns_test(A, test_name);
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
