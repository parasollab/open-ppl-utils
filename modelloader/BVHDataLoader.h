///////////////////////////////////////////////////////////////////////////////////////////
// This file defines a BVH data loader
// the data is broken up into edges per keyframe 


#ifndef _BVHDATLOADER_H_
#define _BVHDATLOADER_H_

#include "ILoadable.h"
//#include <fstream>
#include "InclDefines.h"


class CBVHDataLoader : public ILoadable
{
public:

    //////////////////////////////////////////////////////////////////////////////////////
    CBVHDataLoader() {
       m_GLID = -1;
    }
    void load(ifstream& infile);

    //////////////////////////////////////////////////////////////////////////////////////
    // Implemetation of ILoadable interface
    //////////////////////////////////////////////////////////////////////////////////////
    virtual bool ParseFile(bool silent=false);

    //////////////////////////////////////////////////////////////////////////////////////
    // Implemetation of IModel interface
    //////////////////////////////////////////////////////////////////////////////////////
    const PtVector & GetVertices() const{ return points; }
    const TriVector & GetTriP() const{ return triP; } 
    const TriVector & GetTriN() const{ return triN; }  //triangle normals
    const TriVector & GetTriT() const{ return triT; }  //triangle texture
    const V2Vcetor & GetTextureCoords() const { return textures; }
    const V3Vcetor & GetNormals() const { return normals; }

    PtVector & GetVertices() { return points; }
    TriVector & GetTriP(){ return triP; } 
    TriVector & GetTriN(){ return triN; }  //triangle normals
    TriVector & GetTriT(){ return triT; }  //triangle texture
    V2Vcetor & GetTextureCoords() { return textures; }
    V3Vcetor & GetNormals() { return normals; }

    //EdgeList& GetEdgeList() { return edgeList; }
////////////////////////////////////////////////////////////////////////////////////////////
//
//  Protected Methods and data members
//
////////////////////////////////////////////////////////////////////////////////////////////
protected:
/*
    virtual bool CheckCurrentStatus(bool silent=false);

    bool ReadOBJ(istream& in);
    bool FirstPass(istream& in);
    bool SecondPass(istream& in);
    //CObjGroup& FindGroup(const string& name);
    //CObjGroup& AddGroup(const string& name);
*/

////////////////////////////////////////////////////////////////////////////////////////////
//
//  Private Methods and data members
//
////////////////////////////////////////////////////////////////////////////////////////////
private:

    PtVector points;
    TriVector triP;  //points index
    TriVector triN;  //normal index
    TriVector triT;  //texture index

    V2Vcetor  textures;
    V3Vcetor  normals;
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    //EdgeList edgeList;
};

class CBVHDataInstance {
 public:
   CBVHDataInstance() { }
   string getFileName() { return m_FileName; }
   void setFileName(string fn) { m_FileName = fn; }
   void setInstanceID(int ii) { m_InstanceID = ii; }
   int getInstanceID() { return m_InstanceID; }
   vector<IModel*>& getModels() { return m_Models; }
   void load(string file, double radius, double height); 

   std::string m_FileName;
   vector<IModel*> m_Models;
   int m_InstanceID;
};

class CBVHDataModelFactory {
 public:
   CBVHDataModelFactory() { }
   void addInstance(CBVHDataInstance* bvhd_ins) { m_BVHDataInstances.push_back( bvhd_ins ); }

   int numInstances() { return m_BVHDataInstances.size(); }
   bool hasBVHDataInstance(string filename, int& index);
   CBVHDataInstance* getInstance(int instIndex) { return m_BVHDataInstances[instIndex]; }

   vector<CBVHDataInstance*> m_BVHDataInstances;
};

#endif

