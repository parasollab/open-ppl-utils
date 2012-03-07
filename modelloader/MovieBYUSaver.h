#ifndef _MOVIEBYU_SAVE_H_
#define _MOVIEBYU_SAVE_H_

//#include <string>
//#include <iostream>
//#include <fstream>
#include "InclDefines.h"
//using namespace std;

#include "ModelGraph.h" //for counting edge size
using namespace modelgraph;

#include "IPolygonal.h"

class MoiveBYUSaver {

    typedef IPolygonal::TriVector TriVector;
    typedef IPolygonal::PtVector PtVector;

public:
    bool SaveBYU(const string& filename, const PtVector& pts, const TriVector& tris)
    {
        //open file
        ofstream fout(filename.c_str());
        if( fout.good()==false() ){
            cerr<<"Can't Open file : "<<filename<<endl;
            return false;
        }

        bool r=SaveBYU(fout,pts,tris);
        fout.close();
        return r;
    }

    bool SaveBYU(ostream& out, const PtVector& pts, const TriVector& tris)
    {
        out<<1<<" "<<pts.size()<<" "<<tris.size()<<" "<<tris.size()*3<<"\n";
        out<<1<<" "<<tris.size()<<"\n";
        for( PtVector::const_iterator ip=pts.begin();ip!=pts.end();ip++ )
            out<<(*ip)[0]<<" "<<(*ip)[1]<<" "<<(*ip)[2]<<"\n";
        for( TriVector::const_iterator it=tris.begin();it!=tris.end();it++ )
            out<<(1+(*it)[0])<<" "<<(1+(*it)[1])<<" "<<-(1+(*it)[2])<<"\n";
        return true;
    }

};

#endif //_MOVIEBYU_SAVE_H_

