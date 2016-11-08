#include <fstream>

#include "nonstd/io.h"
#include "nonstd/runtime.h"

using namespace std;

namespace nonstd {

  vector<string>
  read_file(const string& _filename)
  {
    // Open the file.
    ifstream file(_filename);
    assert_msg(file.is_open(),
        "nonstd::read_file error: could not read file '" + _filename + "'!");

    // Read the file.
    string line;
    vector<string> lines;
    while(getline(file, line))
      lines.push_back(line);

    // Close the file and return the contents.
    file.close();
    return lines;
  }

}
