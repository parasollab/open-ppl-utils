#include "glutils/color.h"

#include "nonstd/io.h"
#include "nonstd/numerics.h"

#include <cctype>


namespace glutils {

  /*--------------------------- Construction ---------------------------------*/

  color::
  color(const GLfloat _r, const GLfloat _g, const GLfloat _b, const GLfloat _a)
      : m_rgba{_r, _g, _b, _a}
  { }

  /*----------------- Implicit Conversions to OpenGL Arrays ------------------*/

  color::
  operator GLfloat*() noexcept
  {
    return m_rgba.data();
  }


  color::
  operator const GLfloat*() const noexcept
  {
    return m_rgba.data();
  }


  color::
  operator GLvoid*() noexcept
  {
    return static_cast<GLvoid*>(m_rgba.data());
  }


  color::
  operator const GLvoid*() const noexcept
  {
    return static_cast<const GLvoid*>(m_rgba.data());
  }

  /*------------------------------ Accessors ---------------------------------*/

  GLfloat&
  color::
  operator[](const unsigned short _i) noexcept
  {
    return m_rgba[_i];
  }


  const GLfloat&
  color::
  operator[](const unsigned short _i) const noexcept
  {
    return m_rgba[_i];
  }


  color::iterator
  color::
  begin() noexcept
  {
    return m_rgba.begin();
  }


  color::iterator
  color::
  end() noexcept
  {
    return m_rgba.end();
  }


  color::const_iterator
  color::
  begin() const noexcept
  {
    return m_rgba.begin();
  }


  color::const_iterator
  color::
  end() const noexcept
  {
    return m_rgba.end();
  }

  /*------------------------------ Equality ----------------------------------*/

  bool
  color::
  operator==(const color& _c) const noexcept
  {
    return nonstd::approx(m_rgba[0], _c.m_rgba[0])
        && nonstd::approx(m_rgba[1], _c.m_rgba[1])
        && nonstd::approx(m_rgba[2], _c.m_rgba[2])
        && nonstd::approx(m_rgba[3], _c.m_rgba[3]);
  }


  bool
  color::
  operator!=(const color& _c) const noexcept
  {
    return !(*this == _c);
  }

  /*------------------------------ Weak Ordering -----------------------------*/

  bool
  color::
  operator<(const color& _c) const noexcept
  {
    for(size_t i = 0; i < 4; ++i) {
      if(m_rgba[i] < _c.m_rgba[i])
        return true;
      else if(m_rgba[i] > _c.m_rgba[i])
        return false;
    }
    return false;
  }

  /*----------------------------- Predefined ---------------------------------*/

  const color color::black         = { 0.,  0.,  0.,  1.};
  const color color::dark_grey     = { .1,  .1,  .1,  1.};
  const color color::medium_grey   = { .3,  .3,  .3,  1.};
  const color color::grey          = { .5,  .5,  .5,  1.};
  const color color::light_grey    = { .7,  .7,  .7,  1.};
  const color color::white         = { 1.,  1.,  1.,  1.};

  const color color::brown         = { .3, .15,  0.,  1.};
  const color color::maroon        = {.35,  0., .15,  1.};

  const color color::dark_red      = { .1, .01, .01,  1.};
  const color color::red           = { 1.,  0.,  0.,  1.};

  const color color::orange        = { 1.,  .5,  0.,  1.};

  const color color::yellow        = { 1.,  1.,  0.,  1.};
  const color color::light_yellow  = { 1.,  1.,  .8,  1.};

  const color color::dark_green    = {.01,  .1, .01,  1.};
  const color color::green         = { 0.,  1.,  0.,  1.};
  const color color::goblin_green  = { .1,  .2,  .1,  1.};

  const color color::midnight_blue = {.08, .06,  .2,  1.};
  const color color::blue_grey     = { .1,  .1,  .2,  1.};
  const color color::blue          = { 0.,  0.,  1.,  1.};
  const color color::light_blue    = { .3,  .3,  .7,  1.};
  const color color::cyan          = { .3,  1.,  1.,  1.};

  const color color::magenta       = { 1.,  0.,  1.,  1.};
  const color color::violet        = { .1,  0.,  .1,  1.};

  /*--------------------------------------------------------------------------*/

}

/*---------------------------- ostream overloads -----------------------------*/

std::ostream&
operator<<(std::ostream& _os, const glutils::color& _c)
{
  _os << "{" << _c[0] << ", " << _c[1] << ", " << _c[2] << ", " << _c[3] << "}";
  return _os;
}


std::istream&
operator>>(std::istream& _is, glutils::color& _c)
{
  char buffer;
  bool leadingDelim = false;

  // Read each of the four values.
  for(size_t i = 0; i < 4; ++i)
  {
    // Clear whitespace.
    _is >> std::ws;

    // Check for a delimiter (not a digit or decimal).
    const unsigned char next = _is.peek();
    if(!isdigit(next) and next != u'.')
    {
      _is >> buffer;
      if(i == 0)
        leadingDelim = true;
    }
    _is >> _c[i];
  }

  // If there was a leading delimiter, extract the following delimiter as well.
  if(leadingDelim)
    _is >> std::ws >> buffer;

  return _is;
}

/*---------------------------------- Hasher ----------------------------------*/

size_t
std::hash<glutils::color>::
operator()(const glutils::color& _c) const noexcept
{
  std::hash<GLfloat> hasher;
  return hasher(_c[0] * 1000 + _c[1] * 100 + _c[2] * 10 + _c[3]);
}

/*----------------------------------------------------------------------------*/
