#include "glutils/obj_file.h"

#include <exception>
#include <fstream>

#include "nonstd/string.h"

namespace glutils {

  obj_file&
  obj_file::
  operator>>(triangulated_model& _t)
  {
    std::ifstream file(m_filename);
    if(!file.good())
      throw std::runtime_error("glutils::obj_file::operator>>() error: could "
          "not open the file '" + m_filename + "'.");

    char c1, c2;
    std::string buffer;
    while(file >> c1) {
      file.get(c2);
      if(c2 != ' ')
        continue;
      switch(c1) {
        case 'v': // It's a vertex, Jim.
        {
          vector3f v;
          file >> v;
          _t.add_point(v, true);
          break;
        }
        case 'f': // It's a facet, Jim.
        {
          getline(file, buffer);
          auto facets = nonstd::tokenize(buffer);
          size_t a = std::stoi(nonstd::tokenize(facets[0], "/").front());
          size_t b = std::stoi(nonstd::tokenize(facets[1], "/").front());
          size_t c = std::stoi(nonstd::tokenize(facets[2], "/").front());
          _t.add_facet(a - 1, b - 1, c - 1);
          break;
        }
        default: // It's a comment or unsupported string, Jim.
          getline(file, buffer); // Eat line.
      }
    }

    return *this;
  }


  obj_file&
  obj_file::
  operator<<(const triangulated_model& _t)
  {
    std::ofstream file(m_filename);
    if(!file.is_open())
      throw std::runtime_error("glutils::obj_file::operator<<() error: could "
          "not open the file '" + m_filename + "'.");

    file << m_message.str() << std::endl
         << "# Vertices: " << _t.num_points() << std::endl
         << "# Facets: " << _t.num_facets() << std::endl;

    for(auto p = _t.points_begin(); p != _t.points_end(); ++p)
      file << "v " << (*p)[0] << " " << (*p)[1] << " " << (*p)[2] << std::endl;
    for(auto f = _t.facets_begin(); f != _t.facets_end(); ++f)
      file << "f " << (*f)[0] + 1 << " " << (*f)[1] + 1
           << " "  << (*f)[2] + 1 << std::endl;

    m_message = std::ostringstream();
    return *this;
  }


  std::ostringstream&
  obj_file::
  msg()
  {
    return m_message;
  }

}