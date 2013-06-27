#ifndef ROTATIONCONVERSIONS_H_
#define ROTATIONCONVERSIONS_H_

#include "EulerAngle.h"
#include "Matrix.h"
#include "Quaternion.h"

namespace mathtool {
  //////////////////////////////////////////////////////////////////////////////
  // Conversions to EulerAngles
  //////////////////////////////////////////////////////////////////////////////

  //assignment from a 3x3 rotation matrix.
  inline EulerAngle& convertFromMatrix(EulerAngle& _e, const Matrix3x3& _m) {
    _e.m_beta = atan2(_m[0][2], sqrt(_m[1][2]*_m[1][2] + _m[2][2]*_m[2][2]));
    if(cos(_e.m_beta) > 0) {
      _e.m_alpha = atan2(-_m[1][2], _m[2][2]);
      _e.m_gamma = atan2(-_m[0][1], _m[0][0]);
    }
    else {
      _e.m_alpha = atan2(_m[1][2], -_m[2][2]);
      _e.m_gamma = atan2(_m[0][1], -_m[0][0]);
    }
    return _e;
  }
  //assignment from a quaternion. convert to matrix first then to Euler.
  Matrix3x3& convertFromQuaternion(Matrix3x3&, const Quaternion& _q);
  inline EulerAngle& convertFromQuaternion(EulerAngle& _e, const Quaternion& _q) {
    Matrix3x3 m;
    convertFromQuaternion(m, _q);
    return convertFromMatrix(_e, m);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Conversions to Matrix3x3
  //////////////////////////////////////////////////////////////////////////////

  //assignment from Euler
  inline Matrix3x3& convertFromEuler(Matrix3x3& _m, const EulerAngle& _e) {
    double sa = sin(_e.m_alpha);
    double ca = cos(_e.m_alpha);
    double sb = sin(_e.m_beta);
    double cb = cos(_e.m_beta);
    double sg = sin(_e.m_gamma);
    double cg = cos(_e.m_gamma);
    _m[0][0] = cb*cg;
    _m[0][1] = -cb*sg;
    _m[0][2] = sb;
    _m[1][0] = sa*sb*cg + ca*sg;
    _m[1][1] = -sa*sb*sg + ca*cg;
    _m[1][2] = -sa*cb;
    _m[2][0] = -ca*sb*cg + sa*sg;
    _m[2][1] = ca*sb*sg + sa*cg;
    _m[2][2] = ca*cb;
    return _m;
  }
  //assignment from Quaternion
  inline Matrix3x3& convertFromQuaternion(Matrix3x3& _m, const Quaternion& _q) {
    double  w = _q.m_s;
    double  x = _q.m_v[0];
    double  y = _q.m_v[1];
    double  z = _q.m_v[2];

    _m[0][0] = 1.0 - 2.0*y*y - 2.0*z*z;
    _m[1][0] = 2.0*x*y - 2.0*w*z;
    _m[2][0] = 2.0*x*z + 2.0*w*y;

    _m[0][1] = 2.0*x*y + 2.0*w*z;
    _m[1][1] = 1.0 - 2.0*x*x - 2.0*z*z;
    _m[2][1] = 2.0*y*z - 2.0*w*x;

    _m[0][2] = 2.0*x*z - 2.0*w*y;
    _m[1][2] = 2.0*y*z + 2.0*w*x;
    _m[2][2] = 1.0 - 2.0*x*x - 2.0*y*y;
    return _m;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Conversions to Quaternions
  //////////////////////////////////////////////////////////////////////////////

  //assignment from Euler. First convert to matrix then to quaternion.
  inline Quaternion& convertFromEuler(Quaternion& _q, const EulerAngle& _e) {
    Matrix3x3 m;
    convertFromEuler(m, _e);
    return convertFromMatrix(_q, m);
  }
  //assignment from 3x3 rotation matrix
  inline Quaternion& convertFromMatrix(Quaternion& _q, const Matrix3x3& _m) {
    double  w, x, y, z;

    double wSquare = 0.25*(1.0 + _m[0][0] + _m[1][1] + _m[2][2]);
    if (wSquare > 1e-6){
      w = sqrt(wSquare);
      x = (_m[1][2] - _m[2][1]) / 4.0*w;
      y = (_m[2][0] - _m[0][3]) / 4.0*w;
      z = (_m[0][1] - _m[1][0]) / 4.0*w;
    }
    else {
      w = 0.0;
      double xSquare = -0.5*(_m[1][1] + _m[2][2]);
      if (xSquare > 1e-6){
        x = sqrt(xSquare);
        y = _m[0][1] / 2.0*x;
        z = _m[0][2] / 2.0*x;
      }
      else{
        x = 0.0;
        double ySquare = 0.5*(1.0 - _m[2][2]);
        if (ySquare > 1e-6){
          y = sqrt(ySquare);
          z = _m[1][2] / 2.0*y;
        }
        else{
          y = 0.0;
          z = 1.0;
        }
      }
    }

    _q.m_s = w;
    _q.m_v = Vector3d(x,y,z);
    return _q;
  }

}

#endif
