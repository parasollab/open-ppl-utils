#ifndef NONSTD_IO_H_
#define NONSTD_IO_H_

#include <array>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace nonstd {

  ///@name IO
  ///@{

  /// Read a file line-by-line.
  /// @param[in] _filename The name of the file to read.
  /// @return A vector holding the lines of the file.
  std::vector<std::string> read_file(const std::string& _filename);


  /// Print the contents of a container to a std::string. The container must
  /// support bi-directional iterators.
  /// @param _c The container to print.
  /// @return A string representation of _c's contents.
  template <typename ContainerType>
  std::string
  print_container(const ContainerType& _c)
  {
    std::ostringstream os;
    os << "{";
    auto last = _c.end();
    --last;
    for(auto i = _c.begin(); i != last; ++i)
      os << *i << ", ";
    os << *last << "}";
    return os.str();
  }

  ///@}

}


/*------------------- operator<< overloads for stl containers ----------------*/

namespace std {

  template <typename T1, typename T2>
  ostream&
  operator<<(ostream& _os, const pair<T1, T2>& _p)
  {
    return _os << "(" << _p.first << ", " << _p.second << ")";
  }


  template <typename... T>
  ostream&
  operator<<(ostream& _os, const list<T...>& _l)
  {
    return _os << nonstd::print_container(_l);
  }


  template <typename... T>
  ostream&
  operator<<(ostream& _os, const map<T...>& _m)
  {
    return _os << nonstd::print_container(_m);
  }


  template <typename... T>
  ostream&
  operator<<(ostream& _os, const unordered_map<T...>& _m)
  {
    return _os << nonstd::print_container(_m);
  }


  template <typename... T>
  ostream&
  operator<<(ostream& _os, const set<T...>& _s)
  {
    return _os << nonstd::print_container(_s);
  }


  template <typename... T>
  ostream&
  operator<<(ostream& _os, const unordered_set<T...>& _s)
  {
    return _os << nonstd::print_container(_s);
  }


  template <typename... T>
  ostream&
  operator<<(ostream& _os, const vector<T...>& _v)
  {
    return _os << nonstd::print_container(_v);
  }


  template <typename T, size_t N>
  ostream&
  operator<<(ostream& _os, const array<T, N>& _a)
  {
    return _os << nonstd::print_container(_a);
  }

}

/*----------------------------------------------------------------------------*/

#endif
