#include "MovieBYULoader.h"
#include "ModelFactory.h"
#include "ObjLoader.h"
#include "ModelTool.h"

IModel * CreateModelLoader(const string& file, bool silent)
{
    ILoadable * model=NULL;
    //get extension
    unsigned int pos=file.rfind('.');
    if( pos==string::npos ){
        if(!silent) cerr<<"! Error : Can't Recognize file :"<<file<<endl<<flush;
        return NULL;
    }
    string ext=file.substr(pos+1);
    if( ext=="g" )
        model=new CMovieBYULoader();
    else if( ext=="obj")
        model=new CObjLoader();
    else
        if(!silent) cerr<<"! Error : Can't Recognize extension *."<<ext<<endl<<flush;

    if( model==NULL ) return NULL;
    model->SetDataFileName(file.c_str());
    if( model->ParseFile(silent)==false ){
        delete model;
        return NULL;
    }
    
    return model;
}

vector< IModel* >& CreateModelLoaderVec(CBVHDataModelFactory* modelFactory,
      const string& file, double radius, double height) {
   int instanceIndex=-1;
   if( modelFactory->hasBVHDataInstance( file, instanceIndex ) ) {
      //return pre-loaded vector of models
      //cout << " In CreateModelLoaderVec -- using pre-loaded model vec, instance: " << instanceIndex << endl;
      return modelFactory->getInstance(instanceIndex)->getModels();
   }
   else {
      CBVHDataInstance* bvhDataInstance=new CBVHDataInstance();
      bvhDataInstance->load( string(file), radius, height );
      modelFactory->addInstance( bvhDataInstance );
      return bvhDataInstance->getModels();
   }

   ////////////////////////////////////////////////////////////////////////////
}

IModel * CreateModelLoaderFromPts( vector<Point2d>& boundary ) {

   cout << " --CreateModelLoaderFromPts." << endl;
   //compute pt center -- this is the part that assumes covnex
   Point2d center(0,0);
   for(unsigned int i=0; i<boundary.size(); i++) {
      center[0]+=(boundary[i])[0];
      center[1]+=(boundary[i])[1];
   }
   center[0]/=boundary.size();
   center[1]/=boundary.size();

   //make dummy CMovieBYULoader
   CMovieBYULoader* byuloader = new CMovieBYULoader();
   //add in pts
   typedef vector<Point3d> PtVector;
   PtVector& ptvec = byuloader->GetVertices();
   for(unsigned int i=0; i<boundary.size(); i++) {
      Point3d newPt( (boundary[i])[0], 0, (boundary[i])[1] );
      ptvec.push_back( newPt );
      cout << " pushing pt: " << newPt << endl;
   }
   //add in center
   Point3d center3(center[0],0,center[1]);
   ptvec.push_back( center3 );
   int centerID = ptvec.size()-1;
   ///////////3d//////////////////////////////////////////
   //push the same triangles in for 3d (at predefined height)
   double height =1.0;
   for(unsigned int i=0; i<boundary.size(); i++) {
      Point3d newPt( (boundary[i])[0], height, (boundary[i])[1] );
      ptvec.push_back( newPt );
      //cout << " pushing pt: " << newPt << endl;
   }
   //add in center -- 3d
   Point3d center3_2(center[0],height,center[1]);
   ptvec.push_back( center3_2 );
   int centerID_2 = ptvec.size()-1;
   //cout << " added " << ptvec.size() << " pts "<< endl;

   //make triangles
   typedef Vector<int, 3> Tri;
   typedef vector<Tri>    TriVector;
   TriVector& triP = byuloader->GetTriP();
   for(unsigned int i=1; i<boundary.size(); i++) {
      int id1=i-1;
      int id2=i;
      Tri tri(id1, id2, centerID);
      triP.push_back( tri );
      //cout << "tri size: " << triP.size() << " ids: ["<<id1<<" "<<id2<<" "<<centerID<<"]"<<endl;
      if( i==boundary.size()-1 ) {
	 //push in last one (that wraps around)
	 id1=i;
	 id2=0;
	 Tri tri(id1, id2, centerID);
	 triP.push_back( tri );
	 //cout << "tri size: " << triP.size() << " ids: ["<<id1<<" "<<id2<<" "<<centerID<<"]"<<endl;
      }
   }
   ///////////3d//////////////////////////////////////////
   int index2=boundary.size()+2;
   for(unsigned int i=index2; i<ptvec.size()-1; i++) {
      int id1=i-1;
      int id2=i;
      Tri tri(id1, id2, centerID_2);
      triP.push_back( tri );
      //cout << "tri size: " << triP.size() << " ids: ["<<id1<<" "<<id2<<" "<<centerID_2<<"]"<<endl;
      if( i==ptvec.size()-2 ) {
	 //push in last one (that wraps around)
	 id1=i;
	 id2=boundary.size()+1;
	 Tri tri(id1, id2, centerID_2);
	 triP.push_back( tri );
	 //cout << "tri size: " << triP.size() << " ids: ["<<id1<<" "<<id2<<" "<<centerID_2<<"]"<<endl;
      }
   }

   ComputeFaceNormal( *byuloader );
   //cout << " added " << triP.size() << " triangles "<< endl;
   //return object
   cout << " done!--CreateModelLoaderFromPts." << endl;
   return byuloader;
}
