#include <algorithm>
#include <iostream>
#include <string>

#include "nonstd/vector.h"
#include "nonstd/numerics.h"
#include "nonstd/runtime.h"

using namespace std;
using namespace nonstd;

static const string er = "\ttest_vector error: ";
typedef vector_type<double, 3> vector3;


void test_instantiation() {
  vector3 zero;
  vector3 one{1, 1, 1};
  double t[3] = {2, 2, 2};
  vector3 two = t;
  vector3 a = two;

  assert_msg(zero.size() == 3 &&
      all_of(zero.begin(), zero.end(), [](const double d){return d == 0;}),
      "default construction failed!");
  assert_msg(one.size() == 3 &&
      all_of(one.begin(), one.end(), [](const double d){return d == 1;}),
      "initializer list construction failed!");
  assert_msg(two.size() == 3 &&
      all_of(two.begin(), two.end(), [](const double d){return d == 2;}),
      "array construction failed!");
  assert_msg(a == two, er + "copy constructor failed!");
}

void test_arithmetic() {
  string when = "when testing arithmetic, ";

  vector3 zero,
          one{1, 1, 1},
          a{1, 2, 3},
          b{3, 2, 1},
          four{4, 4, 4},
          i{1, 0, 0},
          j{0, 1, 0},
          k{0, 0, 1};

  assert_msg(zero + one == one, er + when + "identity addition failed!");
  assert_msg(one - zero == one, er + when + "identity subtraction failed!");
  assert_msg(-one == vector3{-1, -1, -1}, er + when + "unary minus failed!");
  assert_msg(a + b == four, er + when + "vector addition failed!");
  assert_msg(four - b == a, er + when + "vector subtraction failed!");
  assert_msg(one * 4 == four, er + when + "scalar multiplication failed!");
  assert_msg(four / 4 == one, er + when + "scalar division failed!");
  assert_msg(approx(a * b, 10.), er + when + "vector dot product failed!");
  assert_msg(approx(i * j, 0.), er + when + "vector dot product failed!");
  assert_msg(i % j == k, er + when + "vector cross product failed!");
}

void test_arithmetic_assignment() {
  string when = "when testing arithmetic assignment, ";

  vector3 zero,
          one{1, 1, 1},
          a{1, 2, 3},
          four{4, 4, 4};

  assert_msg((vector3() += one) == one, er + when + "identity addition failed!");
  assert_msg((vector3{1, 1, 1} -= one) == zero,
      er + when + "identity subtraction failed!");
  assert_msg((vector3{-1, -2, -3} += a) == zero,
      er + when + "vector addition failed!");
  assert_msg((vector3{-1, -2, -3} -= a) == vector3{-2, -4, -6},
      er + when + "vector subtraction failed!");
  assert_msg((vector3{1, 1, 1} *= 4) == four,
      er + when + "scalar multiplication failed!");
  assert_msg((vector3{2, 2, 2} /= 2) == one,
      er + when + "scalar division failed!");
}

void test_magnitudes() {
}

void test_projections() {
}

void test_utilities() {
}

void test_generators() {
  const string when = "when testing generators, ";

  using vec = vector_type<float, 3>;
  vec x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1};
  assert_msg(vec::make_basis(0) == x, er + when + "expected x-basis, but got "
      + to_string(vec::make_basis(0)) + "!");
  assert_msg(vec::make_basis(1) == y, er + when + "expected y-basis, but got "
      + to_string(vec::make_basis(1)) + "!");
  assert_msg(vec::make_basis(2) == z, er + when + "expected z-basis, but got "
      + to_string(vec::make_basis(2)) + "!");
}


int main() {
  test_instantiation();
  test_arithmetic();
  test_arithmetic_assignment();
  test_magnitudes();
  test_projections();
  test_utilities();
  test_generators();

  cerr << "\ttest_vector passed" << endl;
}
