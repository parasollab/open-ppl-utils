#ifndef GLUTILS_GLTRAITS_H_
#define GLUTILS_GLTRAITS_H_

#include <array>

#include <GL/gl.h>

#include "nonstd/vector.h"

namespace glutils {

  ///@name Constants
  ///@{
  /// Precision chosen based on numeric_limits<GLfloat>::digits10.

  static constexpr GLfloat PI    = 3.14159;
  static constexpr GLfloat TWOPI = 6.28319;
  static constexpr GLfloat RadPerDeg = PI / 180;
  static constexpr GLfloat DegPerRad = 180 / PI;

  ///@}
  ///@name Vector Representations
  ///@{

  typedef nonstd::vector_type<GLfloat,  2> vector2f;
  typedef nonstd::vector_type<GLdouble, 2> vector2d;
  typedef nonstd::vector_type<GLfloat,  3> vector3f;
  typedef nonstd::vector_type<GLdouble, 3> vector3d;

  ///@}
  ///@name Transforms
  ///@{

  /// An OpenGL transform matrix. The convention for glutils is opposite of
  /// OpenGL because C/C++ are designed for row-major order. Thus, the transform
  /// matrix is stored in row-major order with the translation elements in the
  /// right-most column.
  typedef std::array<GLfloat, 16> transform;

  /// Apply an OpenGL transform matrix to the current GL stack.
  /// @param[in] _t The transform to apply.
  void apply_transform(const transform& _t) noexcept;

  /// Generate an identity transform.
  transform identity_transform() noexcept;

  ///@}

}

/*-------------------------- Inlined Functions -------------------------------*/

inline
void
glutils::
apply_transform(const glutils::transform& _t) noexcept
{
  glMultTransposeMatrixf(_t.data());
}


inline
glutils::transform
glutils::
identity_transform() noexcept
{
  return glutils::transform{1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1};
}

/*----------------------------------------------------------------------------*/

#endif
