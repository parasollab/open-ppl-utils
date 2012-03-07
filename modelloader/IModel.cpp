#include "IModel.h"
#include "Matrix.h"

void IModel::updateBBX(Point3d& p) {
   //x
   if (p[0]< m_BBX[0]) m_BBX[0] = p[0];
   if (p[0]> m_BBX[1]) m_BBX[1] = p[0];
   //y
   if (p[1]< m_BBX[2]) m_BBX[2] = p[1];
   if (p[1]> m_BBX[3]) m_BBX[3] = p[1];
   //z
   if (p[2]< m_BBX[4]) m_BBX[4] = p[2];
   if (p[2]> m_BBX[5]) m_BBX[5] = p[2];
}
void IModel::resetBBX() {
   m_BBX.clear();
   for(int i=0; i<6; i++) m_BBX.push_back( 1e10 );
   m_BBX[1]*=-1.0;
   m_BBX[3]*=-1.0;
   m_BBX[5]*=-1.0;
}
void IModel::genUniquePts() {
   if( m_UniquePts.size()==0 ) {
      typedef EdgeList::iterator EIT;
      typedef PtVector::iterator PIT;
      EdgeList& el = edgeList;
      for(EIT eit=el.begin(); eit!=el.end(); eit++) {
	 Edge& edge = *eit;
	 Point3d& p1 = edge.first;
	 Point3d& p2 = edge.second;
	 bool addP1=true;
	 bool addP2=true;
	 for(PIT pit=m_UniquePts.begin(); pit!=m_UniquePts.end(); pit++) {
            double dist1 = (p1-(*pit)).norm();
            double dist2 = (p2-(*pit)).norm();
	    if( dist1 < 0.1 ) addP1 = false;
	    if( dist2 < 0.1 ) addP2 = false;
	 }//endfor pit
	 if( addP1 ) m_UniquePts.push_back( p1 );
	 if( addP2 ) m_UniquePts.push_back( p2 );
      }//endfor

   }
}
void IModel::skelBBX(EdgeList& el) { 
   typedef EdgeList::iterator EIT;
   if (m_BBX.size() == 0 ) {
      resetBBX();
      /*
      for(int i=0; i<6; i++) m_BBX.push_back( 1e10 );
      m_BBX[1]*=-1.0;
      m_BBX[3]*=-1.0;
      m_BBX[5]*=-1.0;
      */
   }
   for(EIT eit=el.begin(); eit!=el.end(); eit++) {
      Edge& edge = *eit;
      Point3d& p1 = edge.first;
      Point3d& p2 = edge.second;
      updateBBX(p1);
      updateBBX(p2);
   }//endfor
   m_Center[0] = (m_BBX[0]+m_BBX[1])/2.0;
   m_Center[1] = (m_BBX[2]+m_BBX[3])/2.0;
   m_Center[2] = (m_BBX[4]+m_BBX[5])/2.0;
   cout << " BBX: [" ;
   for(int i=0; i<6; i++) cout<< m_BBX[i] << " ";
   cout << "] center: "<<m_Center<<"\t height: "<<m_BBX[3]-m_BBX[2]<<endl;
   cout << "skel height: " << m_BBX[3]-m_BBX[2] << endl;
}
void IModel::genDir(EdgeList& el1, IModel* m, EdgeList& el2, int index) {
   Edge& first_edge_i = el1[index];
   Point3d p1 = first_edge_i.first;
   Edge& first_edge_j = el2[index];
   Point3d p2 = first_edge_j.first;
   Vector3d dir3d = p2-p1;
   m_Dir = dir3d;
   m_DirectionChange = dir3d;
   m_DistanceChange = m_DirectionChange.norm();
}
double IModel::getBBXRadius() {
   double xrange = fabs(m_BBX[1]-m_BBX[0])/2.0;
   double zrange = fabs(m_BBX[5]-m_BBX[4])/2.0;
   if (xrange > zrange) return xrange; //return smallest
   else return zrange;
}
double IModel::getRadiusToJoint(EdgeList& el, int jointIndex, bool useFirst) {
   //for center: use m_Center
   //find joint pt
   Edge& edge = el[jointIndex];
   Point3d p;
   if( useFirst ) { p = edge.first; }
   else { p = edge.second; }
   
   //radius is computed in xz-plane

   double dist = sqrt( pow(m_Center[0]-p[0],2.0) + pow(m_Center[2]-p[2],2.0) );
   return dist;
}
void IModel::transform(EdgeList& el, double radius, double height) { 
   /*
   typedef EdgeList::iterator EIT;
   Vector3d vecToCenter(m_Center[0],m_Center[1],m_Center[2]);
   Matrix2x2 rotM( cos(m_RotRAD), -1*sin(m_RotRAD), sin(m_RotRAD), cos(m_RotRAD) );
   Vector2d W_dir(0,1);
   
   double bbx_radius = getBBXRadius();
   double bbx_height = fabs(m_BBX[3]-m_BBX[2]);
   double scale = radius/bbx_radius;//model with get scaled by this factor
   cout << " r1: " << radius << " r2: " << bbx_radius << " scale: " << scale << endl;
   m_Dir = scale * m_Dir;
   m_DirectionChange = scale * m_Dir;
   m_DistanceChange *= scale;
   for(EIT eit=el.begin(); eit!=el.end(); eit++) {
      Edge& edge = *eit;
      Point3d& p1 = edge.first;                    //get points
      Point3d& p2 = edge.second;
      Point3d p1_new = p1 + -1.0 * vecToCenter;    //bring back to center
      Point3d p2_new = p2 + -1.0 * vecToCenter;
      Vector2d v1(p1_new[0],p1_new[2]);            //init 2d vec using x,z vals
      Vector2d v2(p2_new[0],p2_new[2]);
      Vector2d v1_new = rotM * v1;                 //rotate
      Vector2d v2_new = rotM * v2;
      Vector3d v1_3d_new( v1_new[0], p1_new[1], v1_new[1] ); //make 3d vectors for scaling
      Vector3d v2_3d_new( v2_new[0], p2_new[1], v2_new[1] ); 
      v1_3d_new = scale * v1_3d_new;               //scale
      v2_3d_new = scale * v2_3d_new;               
      
      p1[0] = v1_3d_new[0];                           //use new vals to update pos
      //p1[1] = v1_3d_new[1];                           //should now be facing z dir 
      p1[2] = v1_3d_new[2];                           
      p2[0] = v2_3d_new[0];              
      //p2[1] = v2_3d_new[1];                                   
      p2[2] = v2_3d_new[2];                                   
   }//endfor
   resetBBX();
   skelBBX(el);
   double y_offset = m_BBX[2];
   if(1)
   for(EIT eit=el.begin(); eit!=el.end(); eit++) {
      Edge& edge = *eit;
      Point3d& p1 = edge.first;                    //get points
      Point3d& p2 = edge.second;
      p1[1] = height* ( p1[1]-y_offset ) / bbx_height;
      p2[1] = height* ( p2[1]-y_offset ) / bbx_height;
   }//endfor
   resetBBX();
   skelBBX(el);
   */
}
