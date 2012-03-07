#ifndef _MODEL_FACTORY_H_
#define _MODEL_FACTORY_H_

//#include <string>
#include "InclDefines.h"

//using namespace std;

#include "BVHDataLoader.h"

class IModel;
IModel * CreateModelLoader(const std::string& file, bool silent=false);

vector<IModel*>& CreateModelLoaderVec(CBVHDataModelFactory* modelFactory, const std::string& file, double radius, double height);

IModel * CreateModelLoaderFromPts( vector<Point2d>& boundary ); 
#endif //_MODEL_FACTORY_H_

