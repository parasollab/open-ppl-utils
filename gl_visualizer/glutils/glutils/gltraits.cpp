#include "glutils/gltraits.h"

namespace glutils {

  void
  apply_transform(const transform& _t) noexcept
  {
    const auto& t = _t.translation();
    const auto& r = _t.rotation();
    GLfloat trans[16] = {r[0][0], r[0][1], r[0][2], t[0],
                         r[1][0], r[1][1], r[1][2], t[1],
                         r[2][0], r[2][1], r[2][2], t[2],
                               0,       0,       0,    1};
    glMultTransposeMatrixf(trans);
  }

}
