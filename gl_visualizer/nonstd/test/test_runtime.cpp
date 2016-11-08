#include <iostream>

#include "nonstd/runtime.h"
#include "nonstd/runtime.h"

using namespace std;
using namespace nonstd;

int main(int _argc, char* _argv[]) {
  // If assert exits the program, the test worked.
  assert_msg(false, "\ttest_runtime passed");
  // Otherwise not so much.
  cerr << "\ttest_runtime failed" << endl;
}
