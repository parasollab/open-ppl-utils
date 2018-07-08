#include <iostream>
#include <sstream>

#include "nonstd/io.h"
#include "nonstd/runtime.h"

using namespace std;
using namespace nonstd;

static const string er = "\ttest_io error: ";


void
test_read_file()
{
  // Read the test file.
  auto lines = read_file("test/test_file.txt");

  string when = "when testing read_file, ";
  vector<string> expected = {"the quick black cat", "jumped over",
      "the lazy dog"};

  // Check number of lines
  assert_msg(lines.size() == 3, er + when + "expected line count = 3, but got "
      "line count = " + to_string(lines.size()));

  // Check each line
  for(size_t i = 0; i < 3; ++i)
    assert_msg(lines[i] == expected[i], er + when + "expected line " +
        to_string(i) + " = '" + expected[i] + "', but got '" + lines[i] + "'!");
}


void
test_std_ostream()
{
  ostringstream out;

  out << pair<int, char>(0, 'a') << endl
      << list<float>{0., 1.1, 2.2} << endl
      << vector<int>{0, 1, 2, 3, 4} << endl
      << array<short, 4>{0, 1, 2, 3} << endl;

  string expected = "(0, a)\n"
                    "{0, 1.1, 2.2}\n"
                    "{0, 1, 2, 3, 4}\n"
                    "{0, 1, 2, 3}\n";

  assert_msg(out.str() == expected, er + "when testing std operator<< helpers, "
      "expected:\n" + expected + "...but got:\n" + out.str());
}


int main() {
  cerr << "\ttesting io..." << flush;

  test_read_file();
  test_std_ostream();

  cerr << "passed" << endl;
}
