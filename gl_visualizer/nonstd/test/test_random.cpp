#include <iostream>
#include <string>

#include "nonstd/random.h"
#include "nonstd/numerics.h"
#include "nonstd/runtime.h"

using namespace std;
using namespace nonstd;

static const string er = "\ttest_random error: ";

void test_drand() {
  // Average 1000 calls to drand and check that its about .5
  double check = 0.;
  for(size_t i = 0; i < 1000; ++i)
    check += drand();
  check /= 1000.;

  assert_msg(approx(check, .5, .03), er + "when testing drand, expected average "
      "of 1000 trials to be ~= .5, but got " + to_string(check));
}


void test_grand() {
  // Average 1000 calls to grand and check that its about .5
  double check = 0.;
  for(size_t i = 0; i < 1000; ++i)
    check += grand();
  check /= 1000.;

  assert_msg(approx(check, .5, .03), er + "when testing grand, expected average "
      "of 1000 trials to be ~= .5, but got " + to_string(check));
}


int main() {
  test_drand();

  cerr << "\ttest_random passed" << endl;
}
