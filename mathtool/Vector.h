#ifndef VECTOR_H_
#define VECTOR_H_

#include "Basic.h"

namespace mathtool{

  //////////////////////////////////////////////////////////////////////////////
  // General D-dimensional Vector
  //////////////////////////////////////////////////////////////////////////////
  template<class T, size_t D>
    class Vector {
      public:

        //construction
        Vector() {
          for(size_t i = 0; i<D; ++i) m_v[i] = T();
        }
        Vector(const Vector& _v) {
          for(size_t i = 0; i<D; ++i) m_v[i] = _v[i];
        }
        Vector(const T _t[D]) {
          for(size_t i = 0; i<D; ++i) m_v[i] = _t[i];
        }

        //assignment
        Vector& operator=(const Vector& _v) {
          for(size_t i = 0; i<D; ++i) m_v[i] = _v[i];
          return *this;
        }
        Vector& operator=(const T _t[D]) {
          for(size_t i = 0; i<D; ++i) m_v[i] = _t[i];
          return *this;
        }

        //access
        typedef T (&arr)[D];
        typedef const T (&carr)[D];
        operator arr() {return m_v;}
        operator carr() const {return m_v;}
        const T* begin() const {return m_v;}
        const T* end() const {return begin()+D;}
        T* begin() {return m_v;}
        T* end() {return begin()+D;}

        //equality
        bool operator==(const Vector& _v) const {
          for(size_t i = 0; i<D; ++i) if(m_v[i] != _v[i]) return false;
          return true;
        }
        //inequality
        bool operator!=(const Vector& _v) const {
          return !(*this == _v);
        }

        //self addition
        Vector& operator+=(const Vector& _v) {
          for(size_t i = 0; i<D; ++i) m_v[i] += _v[i];
          return *this;
        }
        //self subtraction
        Vector& operator-=(const Vector& _v) {
          for(size_t i = 0; i<D; ++i) m_v[i] -= _v[i];
          return *this;
        }
        //self scalar multiply
        Vector& operator*=(const T& _d) {
          for(size_t i = 0; i<D; ++i) m_v[i] *= _d;
          return *this;
        }
        //self scalar divide
        Vector& operator/=(const T& _d) {
          for(size_t i = 0; i<D; ++i) m_v[i] /= _d;
          return *this;
        }
        //self component *
        Vector& operator^=(const Vector& _v) {
          for(size_t i = 0; i<D; ++i) m_v[i] *= _v[i];
          return *this;
        }

        //negation
        Vector operator-() const {
          Vector v;
          for(size_t i = 0; i<D; ++i) v.m_v[i] = -m_v[i];
          return v;
        }
        //addition
        Vector operator+(const Vector& _v) const {
          Vector v(*this);
          return v += _v;
        }
        //subtraction
        Vector operator-(const Vector& _v) const {
          Vector v(*this);
          return v -= _v;
        }
        //scalar multiply
        Vector operator*(const T& _d) const {
          Vector v(*this);
          return v *= _d;
        }
        //scalar divide
        Vector operator/(const T& _d) const {
          Vector v(*this);
          return v /= _d;
        }
        //component *
        Vector operator^(const Vector& _v) const {
          Vector v(*this);
          return v ^= _v;
        }

        //dot product
        T operator*(const Vector& _v) const {
          T dot = 0;
          for(size_t i = 0; i<D; ++i) dot += m_v[i] * _v[i];
          return dot;
        }
        //magnitude
        T norm() const {
          return std::sqrt(normsqr());
        }
        //magnitude squared
        T normsqr() const {
          return (*this)*(*this);
        }
        //normalized vector
        Vector& selfNormalize() {
          T n = norm();
          if(n < std::numeric_limits<T>::epsilon())
            return *this = Vector();
          return *this /= n;
        }
        Vector normalize() const {
          T n = norm();
          if(n < std::numeric_limits<T>::epsilon())
            return Vector();
          return *this / n;
        }

      private:
        T m_v[D];
    };

  //////////////////////////////////////////////////////////////////////////////
  // Specialized 2-dimensional Vector
  //////////////////////////////////////////////////////////////////////////////
  template<class T>
    class Vector<T,2> {
      public:

        //construction
        Vector(const T& _x = T(), const T& _y = T()) {
          m_v[0] = _x; m_v[1] = _y;
        }
        Vector(const Vector& _v) {
          m_v[0] = _v[0]; m_v[1] = _v[1];
        }
        Vector(const T _t[2]) {
          m_v[0] = _t[0]; m_v[1] = _t[1];
        }

        //assignment
        Vector& operator()(const T& _x, const T& _y){
          m_v[0] = _x; m_v[1] = _y;
          return *this;
        }
        Vector& operator=(const Vector& _v) {
          m_v[0] = _v[0]; m_v[1] = _v[1];
          return *this;
        }
        Vector& operator=(const T _t[2]) {
          m_v[0] = _t[0]; m_v[1] = _t[1];
	  return *this;
	}

        //access
        typedef T (&arr)[2];
        typedef const T (&carr)[2];
        operator arr() {return m_v;}
        operator carr() const {return m_v;}
        const T* begin() const {return &m_v[0];}
        const T* end() const {return begin()+2;}
        T* begin() {return &m_v[0];}
        T* end() {return begin()+2;}

        //equality
        bool operator==(const Vector& _v) const {
          return m_v[0] == _v[0] && m_v[1] == _v[1];
        }
        //inequality
        bool operator!=(const Vector& _v) const {
          return !(*this == _v);
        }

        //self addition
        Vector& operator+=(const Vector& _v) {
          m_v[0] += _v[0]; m_v[1] += _v[1];
          return *this;
        }
        //self subtraction
        Vector& operator-=(const Vector& _v) {
          m_v[0] -= _v[0]; m_v[1] -= _v[1];
          return *this;
        }
        //self scalar multiply
        Vector& operator*=(const T& _d) {
          m_v[0] *= _d; m_v[1] *= _d;
          return *this;
        }
        //self scalar divide
        Vector& operator/=(const T& _d) {
          m_v[0] /= _d; m_v[1] /= _d;
          return *this;
        }
        //self component *
        Vector& operator^=(const Vector& _v) {
          m_v[0] *= _v[0]; m_v[1] *= _v[1];
          return *this;
        }

        //negation
        Vector operator-() const {
          return Vector(-m_v[0], -m_v[1]);
        }
        //addition
        Vector operator+(const Vector& _v) const {
          return Vector(m_v[0] + _v[0], m_v[1] + _v[1]);
        }
        //subtraction
        Vector operator-(const Vector& _v) const {
          return Vector(m_v[0] - _v[0], m_v[1] - _v[1]);
        }
        //scalar multiply
        Vector operator*(const T& _d) const {
          return Vector(m_v[0] * _d, m_v[1] * _d);
        }
        //scalar divide
        Vector operator/(const T& _d) const {
          return Vector(m_v[0] / _d, m_v[1] / _d);
        }
        //component *
        Vector operator^(const Vector& _v) const {
          return Vector(m_v[0] * _v[0], m_v[1] * _v[1]);
        }
        //cross product magnitude
        T operator%(const Vector& _v) const {
          return m_v[0]*_v[1] - m_v[1]*_v[0];
        }

        //dot product
        T operator*(const Vector& _v) const {
          return m_v[0]*_v[0] + m_v[1]*_v[1];
        }
        //magnitude
        T norm() const {
          return std::sqrt(normsqr());
        }
        //magnitude squared
        T normsqr() const {
          return (*this)*(*this);
        }
        //normalized vector
        Vector& selfNormalize() {
          T n = norm();
          if(n < std::numeric_limits<T>::epsilon())
            return *this = Vector();
          return *this /= n;
        }
        Vector normalize() const {
          T n = norm();
          if(n < std::numeric_limits<T>::epsilon())
            return Vector();
          return *this / n;
        }
	
	//vector angle
	double getAngle() const {return atan2(m_v[1],m_v[0]);}
	double getAngled() const{return radToRad(getAngle());}
	//rotate vector
	Vector rotate(double _rad) const{ 
	  const double c=cos(_rad), s=sin(_rad);
	  return Vector(m_v[0]*c-m_v[1]*s,m_v[0]*s+m_v[1]*c);
	}
	Vector rotated(double _deg) const {return rotate(degToRad(_deg));}
      
      private:
        T m_v[2];
    };

  //////////////////////////////////////////////////////////////////////////////
  // Specialized 3-dimensional Vector
  //////////////////////////////////////////////////////////////////////////////
  template<class T>
    class Vector<T,3> {
      public:

        //construction
        Vector(const T& _x = T(), const T& _y = T(), const T& _z = T()) {
          m_v[0] = _x; m_v[1] = _y; m_v[2] = _z;
        }
        Vector(const Vector& _v) {
          m_v[0] = _v[0]; m_v[1] = _v[1]; m_v[2] = _v[2];
        }
        Vector(const T _t[3]) {
          m_v[0] = _t[0]; m_v[1] = _t[1]; m_v[2] = _t[2];
        }

        //assignment
        Vector& operator()(const T& _x, const T& _y, const T& _z){
          m_v[0] = _x; m_v[1] = _y; m_v[2] = _z;
          return *this;
        }
        Vector& operator=(const Vector& _v) {
          m_v[0] = _v[0]; m_v[1] = _v[1]; m_v[2] = _v[2];
          return *this;
        }
        Vector& operator=(const T _v[3]) {
          m_v[0] = _v[0]; m_v[1] = _v[1]; m_v[2] = _v[2];
          return *this;
        }

        //access
        typedef T (&arr)[3];
        typedef const T (&carr)[3];
        operator arr() {return m_v;}
        operator carr() const {return m_v;}
        const T* begin() const {return &m_v[0];}
        const T* end() const {return begin()+3;}
        T* begin() {return &m_v[0];}
        T* end() {return begin()+3;}

        //equality
        bool operator==(const Vector& _v) const {
          return m_v[0] == _v[0] && m_v[1] == _v[1] && m_v[2] == _v[2];
        }
        //inequality
        bool operator!=(const Vector& _v) const {
          return !(*this == _v);
        }

        //self addition
        Vector& operator+=(const Vector& _v) {
          m_v[0] += _v[0]; m_v[1] += _v[1]; m_v[2] += _v[2];
          return *this;
        }
        //self subtraction
        Vector& operator-=(const Vector& _v) {
          m_v[0] -= _v[0]; m_v[1] -= _v[1]; m_v[2] -= _v[2];
          return *this;
        }
        //self scalar multiply
        Vector& operator*=(const T& _d) {
          m_v[0] *= _d; m_v[1] *= _d; m_v[2] *= _d;
          return *this;
        }
        //self scalar divide
        Vector& operator/=(const T& _d) {
          m_v[0] /= _d; m_v[1] /= _d; m_v[2] /= _d;
          return *this;
        }
        //self component *
        Vector& operator^=(const Vector& _v) {
          m_v[0] *= _v[0]; m_v[1] *= _v[1]; m_v[2] *= _v[2];
          return *this;
        }
        //self cross product
        Vector& operator%=(const Vector& _v) {
          T v0 = m_v[0], v1 = m_v[1], v2 = m_v[2];
          m_v[0] = v1 * _v[2] - v2 * _v[1];
          m_v[1] = v2 * _v[0] - v0 * _v[2];
          m_v[2] = v0 * _v[1] - v1 * _v[0];
          return *this;
        }

        //negation
        Vector operator-() const {
          return Vector(-m_v[0], -m_v[1], -m_v[2]);
        }
        //addition
        Vector operator+(const Vector& _v) const {
          return Vector(m_v[0] + _v[0], m_v[1] + _v[1], m_v[2] + _v[2]);
        }
        //subtraction
        Vector operator-(const Vector& _v) const {
          return Vector(m_v[0] - _v[0], m_v[1] - _v[1], m_v[2] - _v[2]);
        }
        //scalar multiply
        Vector operator*(const T& _d) const {
          return Vector(m_v[0] * _d, m_v[1] * _d, m_v[2] * _d);
        }
        //scalar divide
        Vector operator/(const T& _d) const {
          return Vector(m_v[0] / _d, m_v[1] / _d, m_v[2] / _d);
        }
        //component *
        Vector operator^(const Vector& _v) const {
          return Vector(m_v[0] * _v[0], m_v[1] * _v[1], m_v[2] * _v[2]);
        }
        //cross product
        Vector operator%(const Vector& _v) const {
          Vector v(*this);
          return v %= _v;
        }

        //dot product
        T operator*(const Vector& _v) const {
          return m_v[0]*_v[0] + m_v[1]*_v[1] + m_v[2]*_v[2];
        }
        //magnitude
        T norm() const {
          return std::sqrt(normsqr());
        }
        //magnitude squared
        T normsqr() const {
          return (*this)*(*this);
        }
        //normalized vector
        Vector& selfNormalize() {
          T n = norm();
          if(n < std::numeric_limits<T>::epsilon())
            return *this = Vector();
          return *this /= n;
        }
        Vector normalize() const {
          T n = norm();
          if(n < std::numeric_limits<T>::epsilon())
            return Vector();
          return *this / n;
        }

      private:
        T m_v[3];
    };

  //////////////////////////////////////////////////////////////////////////////
  // Specialized 4-dimensional Vector
  //////////////////////////////////////////////////////////////////////////////
  template<class T>
    class Vector<T,4> {
      public:

        //construction
        Vector(const T& _x = T(), const T& _y = T(), const T& _z = T(), const T& _w=T()) {
          m_v[0] = _x; m_v[1] = _y; m_v[2] = _z; m_v[3] = _w;
        }
        Vector(const Vector& _v) {
          m_v[0] = _v[0]; m_v[1] = _v[1]; m_v[2] = _v[2]; m_v[3] = _v[3];
        }
        Vector(const T _t[4]) {
          m_v[0] = _t[0]; m_v[1] = _t[1]; m_v[2] = _t[2]; m_v[3] = _t[3];
        }

        //assignment
        Vector& operator()(const T& _x, const T& _y, const T& _z, const T& _w){
          m_v[0] = _x; m_v[1] = _y; m_v[2] = _z; m_v[3] = _w;
          return *this;
        }
        Vector& operator=(const Vector& _v) {
          m_v[0] = _v[0]; m_v[1] = _v[1]; m_v[2] = _v[2]; m_v[3] = _v[3];
          return *this;
        }
        Vector& operator=(const T _t[4]) {
          m_v[0] = _t[0]; m_v[1] = _t[1]; m_v[2] = _t[2]; m_v[3] = _t[3];
          return *this;
        }

        //access
        typedef T (&arr)[4];
        typedef const T (&carr)[4];
        operator arr() {return m_v;}
        operator carr() const {return m_v;}
        const T* begin() const {return &m_v[0];}
        const T* end() const {return begin()+4;}
        T* begin() {return &m_v[0];}
        T* end() {return begin()+4;}

        //equality
        bool operator==(const Vector& _v) const {
          return m_v[0] == _v[0] && m_v[1] == _v[1] && m_v[2] == _v[2] && m_v[3] == _v[3];
        }
        //inequality
        bool operator!=(const Vector& _v) const {
          return !(*this == _v);
        }

        //self addition
        Vector& operator+=(const Vector& _v) {
          m_v[0] += _v[0]; m_v[1] += _v[1]; m_v[2] += _v[2]; m_v[3] += _v[3];
          return *this;
        }
        //self subtraction
        Vector& operator-=(const Vector& _v) {
          m_v[0] -= _v[0]; m_v[1] -= _v[1]; m_v[2] -= _v[2]; m_v[2] -= _v[2];
          return *this;
        }
        //self scalar multiply
        Vector& operator*=(const T& _d) {
          m_v[0] *= _d; m_v[1] *= _d; m_v[2] *= _d; m_v[3] *= _d;
          return *this;
        }
        //self scalar divide
        Vector& operator/=(const T& _d) {
          m_v[0] /= _d; m_v[1] /= _d; m_v[2] /= _d; m_v[3] /= _d;
          return *this;
        }
        //self component *
        Vector& operator^=(const Vector& _v) {
          m_v[0] *= _v[0]; m_v[1] *= _v[1]; m_v[2] *= _v[2]; m_v[3] *= _v[3];
          return *this;
        }
        //negation
        Vector operator-() const {
          return Vector(-m_v[0], -m_v[1], -m_v[2],-m_v[3]);
        }
        //addition
        Vector operator+(const Vector& _v) const {
          return Vector(m_v[0] + _v[0], m_v[1] + _v[1], m_v[2] + _v[2], m_v[3] + _v[3]);
        }
        //subtraction
        Vector operator-(const Vector& _v) const {
          return Vector(m_v[0] - _v[0], m_v[1] - _v[1], m_v[2] - _v[2],  m_v[3] - _v[3]);
        }
        //scalar multiply
        Vector operator*(const T& _d) const {
          return Vector(m_v[0] * _d, m_v[1] * _d, m_v[2] * _d, m_v[3] * _d);
        }
        //scalar divide
        Vector operator/(const T& _d) const {
          return Vector(m_v[0] / _d, m_v[1] / _d, m_v[2] / _d, m_v[3] / _d);
        }
        //component *
        Vector operator^(const Vector& _v) const {
          return Vector(m_v[0] * _v[0], m_v[1] * _v[1], m_v[2] * _v[2], m_v[3] * _v[3]);
        }

        //dot product
        T operator*(const Vector& _v) const {
          return m_v[0]*_v[0] + m_v[1]*_v[1] + m_v[2]*_v[2] + m_v[3]*_v[3];
        }
        //magnitude
        T norm() const {
          return std::sqrt(normsqr());
        }
        //magnitude squared
        T normsqr() const {
          return (*this)*(*this);
        }
        //normalized vector
        Vector& selfNormalize() {
          T n = norm();
          if(n < std::numeric_limits<T>::epsilon())
            return *this = Vector();
          return *this /= n;
        }
        Vector normalize() const {
          T n = norm();
          if(n < std::numeric_limits<T>::epsilon())
            return Vector();
          return *this / n;
        }

      private:
        T m_v[4];
    };

  //////////////////////////////////////////////////////////////////////////////
  // Useful functions. Input/Output and commutativity on multiply
  //////////////////////////////////////////////////////////////////////////////
  //for commutativity of scalar multiply
  template<class T, size_t D>
    inline Vector<T,D> operator*(const T& _d, const Vector<T,D>& _v) {
      return _v*_d;
    }

  template<class T, size_t D>
    inline std::ostream& operator<<(std::ostream& _os, const Vector<T,D>& _v) {
      for(size_t i = 0; i<D; ++i) _os << _v[i] << " ";
      return _os;
    }

  template<class T, size_t D>
    inline std::istream& operator>>(std::istream& _is, Vector<T,D>& _v) {
      for(size_t i=0; i<D; ++i) _is >> _v[i];
      return _is;
    }

  //////////////////////////////////////////////////////////////////////////////
  // Typedef common used vector type
  //////////////////////////////////////////////////////////////////////////////
  typedef Vector<double,2> Vector2d;
  typedef Vector2d Point2d;
  typedef Vector<double,3> Vector3d;
  typedef Vector3d Point3d;
  typedef Vector<double,4> Vector4d;
  typedef Vector4d Point4d;

}

#endif
