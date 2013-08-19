#include "BVHDataLoader.h"
#include "ModelTool.h"
#include <sstream>
#include <cstdio>
//////////////////////////////////////////////////////////////////////////////////////
// Implemetation of ILoadable interface
//////////////////////////////////////////////////////////////////////////////////////
void CBVHDataLoader::load(ifstream& infile) {

   string first;
   int itmp;
   string stmp;

   int keyframenum;
   if( infile >> keyframenum ) { cout << " Reading KF: " << keyframenum << endl; }
   else { cout << " Problem reading keyframe NUM." << endl; exit(-1); }
   while( infile >> first ) {
      if( first == "edge" ) {
	 infile >>itmp;
	 infile>>stmp;
	 Point3d p1,p2;
	 if( stmp=="p1" ) {
	    if(infile>>p1) {
	    }
	    else { cout << "Problem reading p1"<<endl; }
	 }
	 infile>>stmp;
	 if( stmp=="p2" ) {
	    if(infile>>p2) {
	    }
	    else { cout << "Problem reading p2"<<endl; }
	 }
	 Edge e(p1,p2);
	 edgeList.push_back( e );

      }
      else if( first == "endkeyframe" ) {
	 break;
      }
      else {
	 cout << "ERROR: unknown option: " << first <<endl;
      }
   }//endwhile

   setDoSkel(true);

}

bool
CBVHDataLoader::ParseFile(bool silent)
{
   cout << " CBVHDataLoader::ParseFile --> does not ParseFile the way the other loaders do..." << endl;
   return true;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
static int BVH_dataInstanceIndex=0;

void CBVHDataInstance::load(string file, double radius, double height) {
   m_InstanceID=BVH_dataInstanceIndex;//set this instance index
   BVH_dataInstanceIndex++; //increment static index
   setFileName( file );
   // read parameters
   cout << "CBVHDataInstance--> reading: " << file << " instance index: "<<m_InstanceID<< endl;
   char * filename = (char*) file.c_str();
   ifstream infile( filename );
   if( !infile ) {
      cout << "File: " << file << " does not exist. " << endl;
      exit(0);
   }
   string first;
   string stmp;
   //load all models from same file (bvh_data)
   while( infile >> first ) {
      if( first == "keyframe" ) { //load keyframe
	 ILoadable * model=NULL;
	 //get extension
	 unsigned int pos=file.rfind('.');
	 if( pos==string::npos ){
	    cerr<<"! Error : Can't Recognize file :"<<file<<endl;
	    exit(0);
	    //return NULL;
	 }
	 string ext=file.substr(pos+1);
	 model=new CBVHDataLoader();
	 if( model==NULL ) {
	    cout << " mode eq NULL(in bvh_data)...exit!"<<endl;
	    exit(0);
	 }
	 ((CBVHDataLoader*)model)->load(infile);
	 model->skelBBX( ((CBVHDataLoader*)model)->GetEdgeList() );
	 m_Models.push_back( model );
      }//endif
      else {
	 cout << "Should be done reading: " << file << "!!!" << endl;
	 break;
      }
   }//endwhile

   // find transform (might just want a function to find transform
   // need direction, and compute theta
   Vector2d W_dir(0,1);//corresponds to x-z dir, y is for height
   for(int i=0; i<(int)m_Models.size(); i++) {
      IModel* model_i = m_Models[i];
      IModel* model_next;
      if( i==((int)m_Models.size()-1) ) {
	 //just use previous
	 model_i->m_Dir = m_Models[i-1]->m_Dir;
	 Vector2d dir2d(model_i->m_Dir[0], model_i->m_Dir[2]);
	 dir2d=dir2d.normalize();
	 double theta_rad = acos( dir2d*W_dir );
	 model_i->setRotRAD( theta_rad );
	 //double theta_deg = RadToDeg(theta_rad);
	 //cout << "Dir: " << i << "\tdirection: \t" << model_i->m_Dir <<"\t2d dir: "<< dir2d<< endl;
      }
      else {
	 model_next = m_Models[i+1];
	 model_i->genDir(  ((CBVHDataLoader*)model_i)->GetEdgeList(), model_next, ((CBVHDataLoader*)model_next)->GetEdgeList(), 0 );
	 //model_i->m_Dir = model_next->m_Center - model_i->m_Center;
	 Vector2d dir2d(model_i->m_Dir[0], model_i->m_Dir[2]);
	 dir2d=dir2d.normalize();
	 double theta_rad = acos( dir2d*W_dir );
	 model_i->setRotRAD( theta_rad );
	 //double theta_deg = RadToDeg(theta_rad);
	 //cout << "Dir: " << i << "\tdirection: \t" << model_i->m_Dir <<"\t2d dir: "<< dir2d<<" theta: " << theta_deg << endl;
      }
   }//endfor

   // transform
   for(int i=0; i<(int)m_Models.size(); i++) {
      IModel* model_i = m_Models[i];
      model_i->transform(((CBVHDataLoader*)model_i)->GetEdgeList(), radius, height);
   }//endfor

   //return models; //no need to return...Instance has models loaded
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

bool CBVHDataModelFactory::hasBVHDataInstance(string filename, int& index) {

   for(int i=0; i<(int)m_BVHDataInstances.size(); i++) {
      if( m_BVHDataInstances[i]->getFileName() == filename ) {
	 index = i;
	 return true;
      }
   }//endfor i
   index=-1;
   return false;
}
