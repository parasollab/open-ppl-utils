#ifndef _IPOLYGONAL_H_
#define _IPOLYGONAL_H_

//////////////////////////////////////////////////////////////////////
//  This is an interface for all classes which constians
//  polygonal data. For example, deformable object, sculpture deformer
//  and evene model loader
//////////////////////////////////////////////////////////////////////
#include <vector>
using namespace std;

#include <Vector.h>
using namespace mathtool;

class IModel  
{
public:

//////////////////////////////////////////////////////////////////////
//  Interface for Retrive general polgons
//////////////////////////////////////////////////////////////////////
    typedef Vector<int, 3> Tri;

    typedef vector<Point3d>  PtVector;
    typedef vector<Tri>      TriVector;
    typedef vector<Vector3d> V3Vcetor;
    typedef vector<Vector2d> V2Vcetor;
//////////////////////////////////////////////////////////////////////
    virtual const PtVector & GetVertices() const =0;
    virtual const TriVector & GetTriP() const =0;  //triangle points
    virtual const TriVector & GetTriN() const =0;  //triangle normals
    virtual const TriVector & GetTriT() const =0;  //triangle texture

    virtual const V2Vcetor & GetTextureCoords() const =0;
    virtual const V3Vcetor & GetNormals() const =0;
//////////////////////////////////////////////////////////////////////  
    virtual PtVector & GetVertices() =0;
    virtual TriVector & GetTriP() =0;  //triangle points
    virtual TriVector & GetTriN() =0;  //triangle normals
    virtual TriVector & GetTriT() =0;  //triangle texture
    
    virtual V2Vcetor & GetTextureCoords() =0;
    virtual V3Vcetor & GetNormals() =0;

////////////////////////////////////////////////////////////////////// 
//  Skeleton data, typedefs, etc 
////////////////////////////////////////////////////////////////////// 
    typedef pair<Point3d,Point3d> Edge;
    typedef vector<Edge> EdgeList;
    PtVector m_UniquePts;
    EdgeList& GetEdgeList() { return edgeList; }
    void genUniquePts();
    PtVector& GetUniquePts() { return m_UniquePts; }
    EdgeList edgeList;
    bool m_DoSkel;
    void setDoSkel(bool sk) { m_DoSkel = sk; }
    bool getDoSkel() { return m_DoSkel; }
    vector<double> m_BBX;
    double getBBXRadius();
    double getRadiusToJoint(EdgeList& el, int jointIndex, bool useFirst); 
    Point3d m_Center;
    Vector3d m_Dir;
    Vector3d m_DirectionChange;
    double m_DistanceChange;
    double getDistanceChange() { return m_DistanceChange; }
    double m_RotRAD;
    void setRotRAD(double _rotR) { m_RotRAD = _rotR; }
    double getRotRAD() { return m_RotRAD; }
    void resetBBX();
    void updateBBX(Point3d& p);
    void skelBBX(EdgeList& el);
    void genDir(EdgeList& el1, IModel* m, EdgeList& el2, int index);
    void transform(EdgeList& el, double radius, double height);
    int m_GLID;
    int getGLID() { return m_GLID; }
    void setGLID(int glid) { m_GLID=glid; }

};

#endif // !defined(_IPOLYGONAL_H_)

