#ifndef GLUTILS_COLOR_H_
#define GLUTILS_COLOR_H_

#include "glutils/gltraits.h"

namespace glutils {

  //////////////////////////////////////////////////////////////////////////////
  /// Defines various colors in OpenGL's 4fv format. Each color has an red,
  /// green, blue, and alpha value in the range [0,1].
  //////////////////////////////////////////////////////////////////////////////
  class color {

    GLfloat m_rgba[4];   ///< The red, green, blue, alpha values.

    public:

      ///@name Construction
      ///@{

      color(GLfloat _r = 0, GLfloat _g = 0, GLfloat _b = 0, GLfloat _a = 1);

      ///@}
      ///@name Conversion to OpenGL format
      ///@{

      operator const GLfloat*() const noexcept;

      ///@}
      ///@name Accessors
      ///@}

      GLfloat& operator[](const unsigned short _i) noexcept;
      const GLfloat& operator[](const unsigned short _i) const noexcept;

      typedef GLfloat* iterator;
      typedef const GLfloat* const_iterator;

      iterator begin() {return m_rgba;}
      iterator end() {return m_rgba + 4;}
      const_iterator begin() const {return m_rgba;}
      const_iterator end() const {return m_rgba + 4;}

      ///@}
      ///@name Equality
      ///@{

      bool operator==(const color& _c) const noexcept;
      bool operator!=(const color& _c) const noexcept;

      ///@}
      ///@name Predefined Colors
      ///@{

      static const color black;
      static const color dark_grey;
      static const color medium_grey;
      static const color grey;
      static const color light_grey;
      static const color white;

      static const color brown;
      static const color maroon;

      static const color dark_red;
      static const color red;

      static const color orange;

      static const color yellow;
      static const color light_yellow;

      static const color dark_green;
      static const color green;
      static const color goblin_green;

      static const color midnight_blue;
      static const color blue_grey;
      static const color blue;
      static const color light_blue;
      static const color cyan;

      static const color magenta;
      static const color violet;

      ///@}

  };

}

/*---------------------------- ostream overloads -----------------------------*/

std::ostream& operator<<(std::ostream& _os, const glutils::color& _c);

/*----------------------------------------------------------------------------*/

#endif
