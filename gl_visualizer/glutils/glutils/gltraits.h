#ifndef GLUTILS_GLTRAITS_H_
#define GLUTILS_GLTRAITS_H_

#include <GL/gl.h>

#include "nonstd/transform.h"
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

  typedef nonstd::transform_type<GLfloat> transform;

  void apply_transform(const transform& _t) noexcept;

  ///@}

}

#endif
