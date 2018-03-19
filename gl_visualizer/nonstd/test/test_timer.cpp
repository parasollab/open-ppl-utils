#include <iostream>
#include <string>

#include <unistd.h>

#include "nonstd/timer.h"
#include "nonstd/numerics.h"
#include "nonstd/runtime.h"

using namespace std;
using namespace nonstd;

static const string er = "\ttest_timing error: ";


int main() {
  cerr << "\ttesting timing..." << flush;

  timer t;

  assert_msg(approx(t.elapsed(), 0.), er + "expected initial duration ~= 0, but "
      "got " + to_string(t.elapsed()) + "!");

  t.start();
  usleep(22000);
  t.stop();

  assert_msg(approx(t.elapsed(), 22000000., 1000000.), er + "after "
      "usleep(22000), expected duration ~= 2.2e6, but got " +
      to_string(t.elapsed()) + "!");

  cerr << "passed" << endl;
}
