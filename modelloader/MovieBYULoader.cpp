#include <exception>
#include <iostream>
#include <string>
#include <sstream>

#include "MovieBYULoader.h"

using namespace std;


bool
CMovieBYULoader::
ParseFile(bool _silent) {
  if(m_filename == "")
    return false;

  ifstream fin(m_filename);
  if(!fin.good()) {
    if(!_silent)
      cerr << "File " << m_filename << " not found" << endl;
    return false;
  }

  if(!ParseHeader(fin) || !ParseVertices(fin) || !ParsePolygons(fin))
    return false;

  ComputeFaceNormal();
  return true;
}


bool
CMovieBYULoader::
ParseHeader(ifstream& _in) {
  _in >> m_partsSize >> m_vertexSize >> m_polygonSize >> m_edgeSize;
  m_parts.reserve(m_partsSize);

  for(int i = 0; i < m_partsSize && _in; ++i) {
    pair<int,int> range;
    _in >> range.first >> range.second;
    m_parts.push_back(range);
  }
  return true;
}


bool
CMovieBYULoader::
ParseVertices(ifstream& _in) {
  m_points.reserve(m_vertexSize);
  m_cgalPoints.reserve(m_vertexSize);

  string buf;
  getline(_in, buf); // Clear out left-over endline from part specification.

  for(int i = 0; i < m_vertexSize && _in; ++i) {
    getline(_in, buf);
    istringstream is(buf + " " + buf);

    CGALPoint cp;
    Point3d pt;

    is >> cp >> pt;
    m_cgalPoints.push_back(cp);
    m_points.push_back(pt);
  }
  return true;
}


bool
CMovieBYULoader::
ParsePolygons(ifstream& _in) {
  for(int i = 0; i < m_polygonSize; ++i) {
    int id = 0;
    vector<int> index;
    do {
      _in >> id;
      index.push_back(id);
    } while(id >= 0);

    if(index.size() != 3) {
      ostringstream oss;
      oss << "BYU Loader error: when trying to read " << m_filename
          << ", found triangle with " << index.size() << " points on line "
          << 2 + m_partsSize + m_vertexSize + i
          << ". Please trianglate model." << endl;
      throw runtime_error(oss.str());
    }

    index[2] = -index[2];
    Tri tri(index[0] - 1, index[1] - 1, index[2] - 1);
    m_triP.push_back(tri);
  }

  return true;
}
