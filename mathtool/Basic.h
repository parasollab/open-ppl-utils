//----------------------------------------------
//    Definition of general utility routines
//----------------------------------------------
//
//  Programmer: Donald House
//  Date: March 8, 1999
//
//

#ifndef _H_UTILITY
#define _H_UTILITY

#ifndef GB
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>   // define C++ stream I/O routines
#include <iomanip>
using namespace std;
#else
#include "Defines.h"
#endif

//namespace mathtool{

    /* range of real numbers */
    #define SMALLNUMBER 1.0e-10
    #define HUGENUMBER  1.0e10

    /* Miscellaneous Scalar Math */
  
    #define mathtool_abs(x)      (((x) < 0) ? (-(x)) : (x))

    /*
    #define sqr(x)      ((x) * (x))

    #ifndef min
    #define min(x1,x2)  ((x1)<(x2)?(x1):(x2))
    #endif

    #ifndef max
    #define max(x1,x2)  ((x1)>(x2)?(x1):(x2))
    #endif
    
    #define sign(x)     ((x)>=0? 1: -1)
    */
   ///Return minimun between a and b.
   inline double min(double a, double b){
      return a < b ? a : b;
   }

   ///Return maximun between a and b.
   inline double max(double a, double b){
      return a > b ? a : b;
   }

   /// Return the square of a.
   inline double sqr(double a)
   {
      return a*a;
   }

   /// Return sign of x 
   inline int sign(double x) {
      return     ((x)>=0? 1: -1);
   }
      

    //int round(double x, double p);
    //int round( double v );
    inline int round( double x, double p){
        return (int)(((int)((x)*pow(10,p)+((x)<0?-0.5:0.5)))/pow(10,p));
    }

    /*
    inline int round( double v ){
        int integer=(int)floor(v);
        double fraction=v-integer;

        if(v>0)
            return (fraction>=0.5)?integer+1:integer;
        else
            return (fraction>=-0.5)?integer:integer+1;
    }
    */

    //#define swap(x1, x2)  {int tmp=x1; x1=x2; x2=tmp}
    #define applysign(x, y) ((y) >= 0? mathtool_abs(x): -mathtool_abs(x))

    /* Angle Conversions & Constants */

    #ifndef PI
    #define PI 3.1415926535897
    #endif

    #ifndef PI2
    #define PI2 6.2831853071794
    #endif

    #define RAD2DEG (180/PI)
    #define DEG2RAD (PI/180)

    #define DegToRad(x) ((x)*DEG2RAD)
    #define RadToDeg(x) ((x)*RAD2DEG)

    /*
      computes sqrt(a^2 + b^2) without destructive underflow or overflow
    */
    double pythag(double a, double b);

    /*
      Utility Error message routines
    */
    // print s to stdout with trailing blank and no terminating end of line
    void prompt(char *s);

    // print s1, s2, s3 to stdout as blank separated fields with terminating eol
    void message(char *s1, char *s2 = NULL, char *s3 = NULL);

    // print Status: to stdout followed by message(s1, s2, s3)
    void status(char *s1, char *s2 = NULL, char *s3 = NULL);

    // print Error: followed by s1, s2 and s3 to stderr as blank separated fields 
    // with terminating eol
    void error(char *s1, char *s2 = NULL, char *s3 = NULL);

    // print error(s1, s2, s3) and then exit program with code 1 
    void abort(char *s1, char *s2 = NULL, char *s3 = NULL);

    ///Added by Jyh-Ming Lien
    bool getDoubleValue(char * pTag, double * pValue,int size);
    bool getDoubleValue(char * pTag, double * pValue);

    //////////////////////////////////////////////////////////////////////////
    inline void rotateY( double R[3][3], double rY )
    {
        double c=cos(rY);
        double s=sin(rY);
        R[0][0]=c;  R[0][1]=0; R[0][2]=s;
        R[1][0]=0;  R[1][1]=1; R[1][2]=0;
        R[2][0]=-s; R[2][1]=0; R[2][2]=c;
    }


    #ifdef _WIN32

    ////////////////////////////////////////////////////////////////////////////////////////
    // Following functions define M_PI and drand48, which are not starndard c library and 
    // definitions. In addition, rint used to round off float points to int is also here.
    /////////////////////////////////////////////////////////////////////////////////////////

    #define M_PI 3.1415926 //reference PI above

    extern "C" {
        //Implementation of these functions are located in util.cpp
        double drand48();
        double erand48(register unsigned short *xsubi);
        long irand48(register unsigned short m);
        long krand48(register unsigned short *xsubi, unsigned short m);
        long lrand48();
        long mrand48();
        static void next();
        void srand48(long seedval);
        unsigned short * seed48(unsigned short seed16v[3]);
        void lcong48(unsigned short param[7]);
        long nrand48(register unsigned short *xsubi);
        long jrand48(register unsigned short *xsubi);

        /**Round to closest integer.
          *The rint() function rounds x to an integer value according
          *to the prevalent rounding mode.  The default rounding mode
          *is to round to the nearest integer.
          *@return The  rint() function returns the integer value as a float-
          *ing-point number.
          */
        double rint(double x);

    } //end extern "C"

    #endif //_WIN32

//} //end of nprmlib

#endif
