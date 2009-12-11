/********************************************************************

  Vector.H    Header File
  
    Vector Algebra Objects, Methods, and Procedures
    Donald H. House  April 17, 1997
    Visualization Laboratory
    Texas A&M University
    
*********************************************************************/

#ifndef _H_Vector
#define _H_Vector

#include "Basic.h"

#ifndef GB
#include <algorithm>
#include <vector>
using namespace std;
#else
#include "Defines.h"
#endif
////////////////////////////////////////////////////////////////////////////////////////////
//From VectorConstantSize (PMPL)
/**@name 3D rigit body Cfg index*/
//@{
#define Xx  0
#define Yy  1
#define Zz  2
#define ROLLroll    3
#define PITCHpitch  4
#define YAWyaw      5
//@}

//namespace mathtool{
    
    /* Vector Descriptions and Operations */
    template<class T, int D=3 >
    class Vector
    {
    public:

        typedef vector<T> VT;

        Vector(const VT& V);
        Vector(const Vector& V);
        Vector(const T& x=0, const T& y=0, const T& z=0, const T& w=0);
	Vector(std::istream& );
	//next two may be needed for compatibility
        //Vector3(const T& x=0, const T& y=0, const T& z=0);
        Vector(const T& x, const T& y, const T& z, const T& a, const T& b, const T& c);

        T& operator[](int i){ return v[i]; }
        const T& operator[](int i) const{ return v[i]; }
        
        T norm() const;               // magnitude of vector
        T normsqr() const;            // magnitude squared
        Vector normalize() const;     // normalize
	//void selfNormalize();         // normalize this vector
        
        void set(const Vector &other){ *this=other; }
        void set(const VT& V){ v=V; }
        void set(const T other[D]){ 
            /*copy(&other[0],&other[4],v.begin());*/ 
            for( int i=0;i<D;i++ ) v[i]=other[i];
        }
        void set(const T& x=0, const T& y=0, const T& z=0, const T& w=0){
            if( D>0 ) v[0]=x; if( D>1 ) v[1]=y;
            if( D>2 ) v[2]=z; if( D>3 ) v[3]=w;
        }
        void get(T other[D]) const { copy(v.begin(),v.end(),other); }
        void get(VT & other) const { other=v; }

        /* Vector operator prototypes */
        const Vector& operator=(const Vector& v2);        // assignment
        Vector operator-() const;                 // unary negation
        Vector operator+(const Vector& v2) const; // vector add
        Vector operator-(const Vector& v2) const; // vector sub
        Vector operator*(const T& s) const;       // scalar multiply
        Vector operator/(const T& s) const;       // division by scalar
        
        Vector operator^(const Vector& v2) const; // component *
        T operator*(const Vector& v2) const;      // dot product
        Vector operator%(const Vector& v2) const; // cross product
        bool operator==(const Vector& two) const; // equality

	///////////////////////////////////////////////////////////////////////////////
	//functions from PMPL VectorConstantSize
	/// Dot Product Operation. Vnew = (V11*V21)+(V12*V22)+(V13*V23)+....
	//inline double dotProduct(Vector&) const;
	T dotProduct(Vector&) const;
	T magnitude();    
	/// I/O for vector
	void Read(std::istream&);
	void Write(std::ostream&) const;
	Vector crossProduct(const Vector& v2) const; //call %
	T getX() const; 
	T getY() const;
	T getZ() const;
	T getRoll() const;
	T getPitch() const;
	T getYaw() const;
        
    protected:
        VT v;
    };
    
    //Implementation
    template<class T, int D >Vector<T,D>::
    Vector(const VT& V):v(D,0)
    {
        set(V);
    }

    template<class T, int D >Vector<T,D>::
    Vector(const T& x, const T& y, const T& z, const T& w)
    :v(D,0)
    {
        set(x,y,z,w);
    }
    
    template<class T, int D >Vector<T,D>::
    Vector(const T& x, const T& y, const T& z, const T& a, const T& b, const T& c)
    :v(D,0)
    {
            if( D>0 ) v[0]=x; if( D>1 ) v[1]=y;
            if( D>2 ) v[2]=z; if( D>3 ) v[3]=a;
            if( D>4 ) v[4]=b; if( D>5 ) v[5]=c;
    }

    template<class T, int D >
    Vector<T,D>::Vector(const Vector<T,D>& V)
    :v(D,0)
    {
        v=V.v;
    }
    
    template<class T, int D >
    Vector<T,D>::
    Vector(std::istream& is)
    :v(D,0)
    {
       Read(is); 
    }
    
    template<class T, int D >
    T Vector<T,D>::norm() const               // magnitude of vector
    {
        return (T)sqrt(normsqr());
    }
    
    template<class T, int D >
    T Vector<T,D>::normsqr() const            // magnitude squared
    {
        T sumsqr = 0;
        for(int i = 0; i < D; i++)
            sumsqr += pow(v[i],2);
        return sumsqr;
    }
    
    template<class T, int D >
    Vector<T,D> Vector<T,D>::normalize() const     // normalize
    {
        Vector<T,D> newv;
        double magnitude = norm();
        int i;
        for(i=0; i<D; i++){
            newv.v[i] = v[i] / magnitude;
            if(fabs(v[i]) > magnitude * HUGENUMBER){
                cerr << "taking the norm of a zero" << D << " Vector" << endl;
                break;
            }
        }
        for(; i < D; i++){
            newv.v[i] = v[i] / magnitude;
        }
        return newv;
    }
   
    /* 
    template<class T, int D >
    void Vector<T,D>::selfNormalize()           // normalize this vector
    {
        double magnitude = norm();
        for(int i=0; i<D; i++)
            v[i] = v[i] / magnitude;
    }
    */
    
    template<class T, int D >
    const Vector<T,D>& Vector<T,D>::operator=(const Vector<T,D>& v2)
    {
        v=v2.v;
        return *this;
    }
    
    template<class T, int D >
    Vector<T,D> Vector<T,D>::operator-() const
    {
        Vector<T,D> newv(*this);
        for( int i=0;i<D;i++ ) newv.v[i]=-newv.v[i];
        return newv;
    }
    
    template<class T, int D >
    Vector<T,D> Vector<T,D>::operator+(const Vector& v2) const
    {
        Vector<T,D> newv(*this);
        for( int i=0;i<D;i++ ) newv.v[i]+=v2.v[i];
        return newv;
    }
    
    template<class T, int D >
    Vector<T,D> Vector<T,D>::operator-(const Vector& v2) const
    {
        Vector<T,D> newv(*this);
        for( int i=0;i<D;i++ ) newv.v[i]-=v2.v[i];
        return newv;
    }
    
    template<class T, int D >
    Vector<T,D> Vector<T,D>::operator*(const T& s) const
    {
        Vector<T,D> newv(*this);
        for( int i=0;i<D;i++ ) newv.v[i]*=s;
        return newv;
    }
    
    template<class T, int D >
    Vector<T,D> Vector<T,D>::operator/(const T& s) const
    {
        Vector<T,D> newv(*this);
        for( int i=0;i<D;i++ ) newv.v[i]/=s;
        return newv;
    }
    
    template<class T, int D >
    Vector<T,D> Vector<T,D>::operator^(const Vector<T,D>& v2) const
    {
        Vector<T,D> newv(*this);
        for( int i=0;i<D;i++ ) newv.v[i]*=v2.v[i];
        return newv;
    }
    
    template<class T, int D >
    T Vector<T,D>::operator*(const Vector<T,D>& v2) const
    {
        T dot=0;
        for( int i=0;i<D;i++ ) dot+=(v[i]*v2.v[i]);
        return dot;
    }
    
    template<class T, int D >
    Vector<T,D> Vector<T,D>::operator%(const Vector<T,D>& v2) const
    {
        if( D>3 ){
            cerr << "cannot take cross product of " << D << "D Vector";
            exit(1);
        }
        if( D==2 ) return Vector<T,D>();
        //D==3
        Vector<T,D> newv(*this);
        newv.v[0]=v[1] * v2.v[2] - v[2] * v2.v[1];
        newv.v[1]=v[2] * v2.v[0] - v[0] * v2.v[2];
        newv.v[2]=v[0] * v2.v[1] - v[1] * v2.v[0];
        return newv;
    }
    
    template<class T, int D >
    bool Vector<T,D>::operator==(const Vector<T,D>& other) const
    {
        return v==other.v;
    }

    template<class T, int D >
    Vector<T,D> operator*(const T& s, const Vector<T,D>& v)
    {
        Vector<T,D> newv(v);
        for( int i=0;i<D;i++ ) newv[i]=newv[i]*s;
        return newv;
    }

    template<class T, int D>
    ostream & operator<<(ostream & out, const Vector<T,D> & v) {
        for( int d=0;d<D;d++ ) out<<v[d]<<" ";
        return out;
    }

    template<class T, int D>
    istream & operator>>(istream & in, Vector<T,D> & vec) {
        T v[D];
        for(int i=0;i<D;i++) in>>v[i];
        vec.set(v);
        return in;
    }

    /////////////////////////////////////////////////////////////////
    // Function from PMPL VectorConstantSize/Vectors
    template<class T, int D >
    T Vector<T,D>::dotProduct(Vector<T,D>& v2) const
    {
        T dot=0;
	dot = (*this)*v2;
        return dot;
    }
    template<class T, int D >
    T Vector<T,D>::magnitude() 
    {
        return norm();
    }
    template<class T, int D >
    void Vector<T,D>::Read(std::istream & _is)
    {
       _is>>(*this);
    }
    template<class T, int D >
    void Vector<T,D>::Write(std::ostream & _os)const
    {
       _os<<(*this);
    }
    template<class T, int D >
    Vector<T,D> Vector<T,D>::crossProduct(const Vector<T,D>& v2) const
    {
       return (*this)%v2;
    }
    template<class T, int D >
    T Vector<T,D>::getX() const
    {
       return v[Xx];
    }
    template<class T, int D >
    T Vector<T,D>::getY() const
    {
       return v[Yy];
    }
    template<class T, int D >
    T Vector<T,D>::getZ() const
    {
       return v[Zz];
    }
    template<class T, int D >
    T Vector<T,D>::getRoll() const
    {
       return v[ROLLroll];
    }
    template<class T, int D >
    T Vector<T,D>::getPitch() const
    {
       return v[PITCHpitch];
    }
    template<class T, int D >
    T Vector<T,D>::getYaw() const
    {
       return v[YAWyaw];
    }
   

    
    /* Typedef common used vector type */
    typedef Vector<double,2> Vector2d;
    typedef Vector<double,3> Vector3d;
    typedef Vector<double,4> Vector4d;

    //for PMPL code
    typedef Vector3d Vector3D;
    typedef Vector<double,6> Vector6D;

//} //end of nprmlib namespace

#endif
