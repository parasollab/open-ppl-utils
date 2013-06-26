/* Defines 3D orientations through rotation matrices. Handles conversion between
 * types of rotations. Input and output is in Euler Angles. 
 * 
 * Many kinds of representation are implemented here, such as
 * - Euler Angle
 * - Matrix
 * - Quaternion
 * Operations for their representations are provided.
 * Coversion between these representation types are also available.
 */

#ifndef ORIENTATION_H_
#define ORIENTATION_H_

#include "RotationConversions.h"
#include "Vector.h"

namespace mathtool {

  class Orientation {
    public:

      // The type of Orientation instance. 
      enum OrientationType {
        EULER,
        MATRIX,
        QUATERNION
      };

      //Constructor for EulerAngles
      //Also default constructor
      Orientation(const EulerAngle& _e = EulerAngle()) : m_type(EULER), m_euler(_e) {}

      //Constructor for Matrix
      Orientation(const Matrix3x3& _m) : m_type(MATRIX), m_matrix(_m) {}

      //Constructor for Quaternions
      Orientation(const Quaternion& _q) : m_type(QUATERNION), m_quaternion(_q) {}

      //equality - if they aren't the same type: convert _o to type of *this and
      //compare
      bool operator==(const Orientation& _o) const {
        Orientation o(_o);
        o.convertToType(m_type);
        switch(m_type) {
          case(EULER) :
            return m_euler == o.m_euler;
          case(MATRIX) :
            return m_matrix == o.m_matrix;
          case(QUATERNION) :
            return m_quaternion == o.m_quaternion;
        }
      }
      //inequality
      bool operator!=(const Orientation& _o) const {
        return !(*this == _o);
      }

      //Matrix times vector. If type is Euler convert to matrix then multiply.
      //Both quaternion and matrix apply directly the rotation.
      Vector3d operator*(const Vector3d & _v) {
        switch(m_type) {
          case(EULER) :
            convertToType(MATRIX);
          case(MATRIX) :
            return m_matrix * _v;
          case(QUATERNION) :
            m_quaternion.normalize();
            return (m_quaternion * _v * -m_quaternion).imaginary();
        }
      }

      //Self Multiplication.
      //convert the orientations to the same type and perform composition.
      Orientation& operator*=(const Orientation& _o) {
        Orientation o(_o);
        o.convertToType(m_type);
        switch(m_type) {
          case(EULER) :
            m_euler += o.m_euler; break;
          case(MATRIX) :
            m_matrix = m_matrix * o.m_matrix; break;
          case(QUATERNION) :
            m_quaternion *= o.m_quaternion; break;
        }
        return *this;
      }
      //Orientations can be added together if they are EulerAngles. 
      //convert to Euler angles then perform addition.
      Orientation& operator+=(Orientation& _o) {
        convertToType(EULER);
        _o.convertToType(EULER);
        m_euler += _o.m_euler;
        return *this;
      }
      //Orientations can be subtracted if they are EulerAngles. 
      //convert to Euler angles then perform addition.
      Orientation& operator-=(Orientation& _o) {
        convertToType(EULER);
        _o.convertToType(EULER);
        m_euler -= _o.m_euler;
        return *this;
      }

      //inverse
      Orientation operator-() const {
        Orientation o;
        switch(m_type) {
          case(EULER) :
            o.m_euler = -m_euler;
          case(MATRIX) :
            o.m_matrix = m_matrix.transpose();
          case(QUATERNION) :
            o.m_quaternion = -m_quaternion;
        }
        return o;
      }
      //multiplication
      Orientation operator*(const Orientation& _o) const {
        Orientation o(*this);
        o *= _o;
        return o;
      }
      //addition
      Orientation operator+(Orientation& _o) const {
        Orientation o(*this);
        o += _o;
        return o;
      }
      //subtraction
      Orientation operator-(Orientation& _o) const {
        Orientation o(*this);
        o -= _o;
        return o;
      }

      friend std::istream& operator>>(std::istream& _is, Orientation& _o);
      friend std::ostream& operator<<(std::ostream& _os, const Orientation& _o);

    private:
      /** convert the Orientation instance from current type to specified type.
       * _newType The type this instance will be converted to.
       */
      void convertToType(OrientationType _newType) {
        if(m_type == _newType) return;
        switch(m_type) {
          case(EULER) :
            switch(_newType) {
              case(MATRIX) :
                convertFromEuler(m_matrix, m_euler);
                return;
              case(QUATERNION) :
                convertFromEuler(m_quaternion, m_euler);
            }
          case(MATRIX) :
            switch(_newType) {
              case(EULER) :
                convertFromMatrix(m_euler, m_matrix);
                return;
              case(QUATERNION) :
                convertFromMatrix(m_quaternion, m_matrix);
            }
          case(QUATERNION) :
            switch(_newType) {
              case(EULER) :
                convertFromQuaternion(m_euler, m_quaternion);
                return;
              case(MATRIX) :
                convertFromQuaternion(m_matrix, m_quaternion);
            }
        }
      }

      OrientationType m_type;  //type of orientation instance
      EulerAngle m_euler;      //storage if type is EULER
      Matrix3x3 m_matrix;      //storage if type is MATRIX
      Quaternion m_quaternion; //storage if type is QUATERNION
  };

  inline std::istream& operator>>(std::istream& _is, Orientation& _o) {
    _o.m_type == Orientation::EULER;
    return _is >> _o.m_euler;
  }

  inline std::ostream& operator<<(std::ostream& _os, const Orientation& _o) {
    Orientation o(_o);
    o.convertToType(Orientation::EULER);
    return _os << _o.m_euler;
  }

}

#endif

