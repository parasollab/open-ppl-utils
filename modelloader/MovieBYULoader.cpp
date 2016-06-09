#include <iostream>

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

  if(ParseSection1(fin) == false)
    return false;
  if(ParseSection2(fin) == false)
    return false;
  if(ParseSection3(fin) == false)
    return false;
  fin.close();

  ComputeFaceNormal();
  return true;
}


bool
CMovieBYULoader::
ParseSection1(ifstream& _in) {
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
ParseSection2(ifstream& _in) {
  m_points.reserve(m_vertexSize);

  for(int i = 0; i < m_vertexSize && _in; ++i) {
    Point3d pt;
    _in >> pt;
    m_points.push_back(pt);
  }
  return true;
}


bool
CMovieBYULoader::
ParseSection3(ifstream& _in) {
  for(int i = 0; i < m_polygonSize; ++i) {
    int id = 0;
    vector<int> index;
    do {
      _in >> id;
      index.push_back(id);
    } while(id >= 0);

    if(index.size() != 3) {
      cerr << "Please trianglate model first" << endl;
      exit(1);
    }

    index[2] = -index[2];
    Tri tri(index[0] - 1, index[1] - 1, index[2] - 1);
    m_triP.push_back(tri);
  }

  return true;
}
