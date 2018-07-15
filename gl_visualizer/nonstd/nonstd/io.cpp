#include <fstream>

#include "nonstd/io.h"
#include "nonstd/runtime.h"

using namespace std;

namespace nonstd {

  vector<string>
  read_file(
      const string& _filename
  ) {
    // Open the file.
    ifstream file(_filename);
    if(!file.is_open())
      throw nonstd::exception(WHERE) << "Could not read file '" << _filename
                                     << "'.";

    // Read the file.
    string line;
    vector<string> lines;
    while(getline(file, line))
      lines.push_back(line);

    // Close the file and return the contents.
    file.close();
    return lines;
  }


  char
  matching_bracket(
      const char _bracket
  ) {
    switch(_bracket)
    {
      case '(':
        return ')';
      case ')':
        return '(';
      case '[':
        return ']';
      case ']':
        return '[';
      case '{':
        return '}';
      case '}':
        return '{';
      case '<':
        return '>';
      case '>':
        return '<';
      default:
        throw nonstd::exception(WHERE) << "Unrecognized bracket '" << _bracket
                                       << "', choices are (){}[]<>.";
    }
  }

}
