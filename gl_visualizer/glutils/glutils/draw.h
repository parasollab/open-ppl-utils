#ifndef GLUTILS_DRAW_H_
#define GLUTILS_DRAW_H_

#include <cmath>

#include "glutils/gltraits.h"


namespace glutils {

  //////////////////////////////////////////////////////////////////////////////
  /// Simple routines for drawing primitive shapes in OpenGL.
  //////////////////////////////////////////////////////////////////////////////
  namespace draw {

    ///@name 2D Primitives
    ///@{

    ////////////////////////////////////////////////////////////////////////////
    /// Draw a pseudo-circle with _segments sides on the x-y plane.
    inline void
    circle(const GLfloat _radius, const unsigned short _segments = 20)
    {
      GLfloat incr = 2. * PI / _segments;
      GLfloat x, y;

      glBegin(GL_TRIANGLE_FAN);
      glVertex2f(0., 0.);
      for(short i = 0; i <= _segments; ++i) {
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex2f(x, y);
      }
      glEnd();
    }


    ////////////////////////////////////////////////////////////////////////////
    /// Draw a wire pseudo-circle.
    inline void
    circle_frame(const GLfloat _radius, const unsigned short _segments = 20)
    {
      GLfloat incr = 2 * PI / _segments;
      GLfloat x, y;

      glBegin(GL_LINE_LOOP);
      for(short i = 0; i <= _segments; ++i) {
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex2f(x, y);
      }
      glEnd();
    }

    ///@}
    ///@name 3D Primitives
    ///@{

    ////////////////////////////////////////////////////////////////////////////
    /// Draw a sphere centered at the origin with axis along the z direction.
    inline void
    sphere(const GLfloat _radius, const unsigned short _segments = 20)
    {
      GLfloat oIncr = 2 * PI / _segments; // Angle increment for x,y coords.
      GLfloat zIncr = PI / _segments;     // Angle increment for z coords.
      GLfloat x, y, z, r;

      // Draw +zHat cap.
      glBegin(GL_TRIANGLE_FAN);
      glVertex3f(0, 0, _radius); // The +zHat pole.
      z = _radius * cos(zIncr);
      r = _radius * sin(zIncr);
      for(short i = 0; i <= _segments; ++i) {
        x = r * cos(oIncr * i);
        y = r * sin(oIncr * i);
        glVertex3f(x, y, z);
      }
      glEnd();

      // Draw main surface.
      GLfloat z2, r2;
      for(short j = 1; j < _segments; ++j) {
        glBegin(GL_TRIANGLE_STRIP);
        z  = _radius * cos(zIncr * j);
        r  = _radius * sin(zIncr * j);
        z2 = _radius * cos(zIncr * (j + 1));
        r2 = _radius * sin(zIncr * (j + 1));
        for(short i = 0; i <= _segments; ++i) {
          x = cos(oIncr * i);
          y = sin(oIncr * i);
          glVertex3f(x * r , y * r ,  z);
          glVertex3f(x * r2, y * r2, z2);
        }
        glEnd();
      }

      // Draw -zHat cap.
      glBegin(GL_TRIANGLE_FAN);
      glVertex3f(0, 0, -_radius);
      z = _radius * cos(zIncr * (_segments - 1));
      r = _radius * sin(zIncr * (_segments - 1));
      for(short i = _segments; i >= 0; --i) {
        x = r * cos(oIncr * i);
        y = r * sin(oIncr * i);
        glVertex3f(x, y, z);
      }
      glEnd();
    }


    ////////////////////////////////////////////////////////////////////////////
    /// Draw a wire sphere centered at the origin with axis along the z direction.
    inline void
    sphere_frame(const GLfloat _radius, const unsigned short _segments = 20)
    {
      GLfloat oIncr = 2 * PI / _segments; // Angle increment for x,y coords.
      GLfloat zIncr = PI / _segments;     // Angle increment for z coords.
      GLfloat x, y, z, r;

      // Draw latitude lines.
      glBegin(GL_LINES);

      // Draw +zHat cap.
      z = _radius * cos(zIncr);
      r = _radius * sin(zIncr);
      for(short i = 0; i < _segments; ++i) {
        glVertex3f(0, 0, _radius);
        x = r * cos(oIncr * i);
        y = r * sin(oIncr * i);
        glVertex3f(x, y, z);
      }

      // Draw main surface.
      GLfloat z2, r2;
      for(short j = 1; j < _segments; ++j) {
        z  = _radius * cos(zIncr * j);
        r  = _radius * sin(zIncr * j);
        z2 = _radius * cos(zIncr * (j + 1));
        r2 = _radius * sin(zIncr * (j + 1));
        for(short i = 0; i <= _segments; ++i) {
          x = cos(oIncr * i);
          y = sin(oIncr * i);
          glVertex3f(x * r , y * r , z);
          glVertex3f(x * r2, y * r2, z2);
        }
      }

      // Draw -zHat cap.
      z = _radius * cos(zIncr * (_segments - 1));
      r = _radius * sin(zIncr * (_segments - 1));
      for(short i = _segments; i > 0; --i) {
        x = r * cos(oIncr * i);
        y = r * sin(oIncr * i);
        glVertex3f(x, y, z);
        glVertex3f(0, 0, -_radius);
      }
      glEnd();

      // Draw longitude lines.
      for(short i = 1; i < _segments; ++i) {
        glPushMatrix();
        z = _radius * cos(zIncr * i);
        r = _radius * sin(zIncr * i);

        glTranslatef(0, 0, z);
        circle_frame(r, _segments);
        glPopMatrix();
      }
    }


    ////////////////////////////////////////////////////////////////////////////
    /// Draw a wire box.
    inline void
    box_frame(const GLfloat _lenX, const GLfloat _lenY, const GLfloat _lenZ)
    {
      // Get half-lengths.
      GLfloat hlX = _lenX / 2.;
      GLfloat hlY = _lenY / 2.;
      GLfloat hlZ = _lenZ / 2.;

      // Draw.
      glBegin(GL_LINE_STRIP);
      glVertex3f( hlX,  hlY,  hlZ);
      glVertex3f(-hlX,  hlY,  hlZ);
      glVertex3f(-hlX,  hlY, -hlZ);
      glVertex3f( hlX,  hlY, -hlZ);
      glVertex3f( hlX,  hlY,  hlZ);
      glVertex3f( hlX, -hlY,  hlZ);
      glVertex3f(-hlX, -hlY,  hlZ);
      glVertex3f(-hlX, -hlY, -hlZ);
      glVertex3f( hlX, -hlY, -hlZ);
      glVertex3f( hlX, -hlY,  hlZ);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(-hlX,  hlY,  hlZ);
      glVertex3f(-hlX, -hlY,  hlZ);
      glVertex3f(-hlX,  hlY, -hlZ);
      glVertex3f(-hlX, -hlY, -hlZ);
      glVertex3f( hlX,  hlY, -hlZ);
      glVertex3f( hlX, -hlY, -hlZ);
      glEnd();
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Draw a box.
    inline void
    box(const GLfloat _lenX, const GLfloat _lenY, const GLfloat _lenZ)
    {
      // Get half-lengths.
      GLfloat hlX = _lenX / 2.;
      GLfloat hlY = _lenY / 2.;
      GLfloat hlZ = _lenZ / 2.;

      // Draw.
      glBegin(GL_TRIANGLE_FAN);
      glVertex3f( hlX,  hlY,  hlZ);
      glVertex3f( hlX, -hlY,  hlZ);
      glVertex3f( hlX, -hlY, -hlZ);
      glVertex3f( hlX,  hlY, -hlZ);
      glVertex3f(-hlX,  hlY, -hlZ);
      glVertex3f(-hlX,  hlY,  hlZ);
      glVertex3f(-hlX, -hlY,  hlZ);
      glVertex3f( hlX, -hlY,  hlZ);
      glEnd();
      glBegin(GL_TRIANGLE_FAN);
      glVertex3f(-hlX, -hlY, -hlZ);
      glVertex3f(-hlX,  hlY, -hlZ);
      glVertex3f( hlX,  hlY, -hlZ);
      glVertex3f( hlX, -hlY, -hlZ);
      glVertex3f( hlX, -hlY,  hlZ);
      glVertex3f(-hlX, -hlY,  hlZ);
      glVertex3f(-hlX,  hlY,  hlZ);
      glVertex3f(-hlX,  hlY, -hlZ);
      glEnd();
    }


    ////////////////////////////////////////////////////////////////////////////
    /// Draw a cone with the base centered at the origin and tip pointed away
    /// from the camera.
    /// pointing at the camera.
    /// @param[in] _radius The radius of the base.
    /// @param[in] _height The height of the cone.
    /// @param[in] _segments The number of segments to use for the sides.
    inline void
    cone(const GLfloat _radius, const GLfloat _height,
        const unsigned short _segments = 20)
    {
      circle(_radius, _segments);

      GLfloat incr = 2 * PI / _segments;
      GLfloat x, y;

      glBegin(GL_TRIANGLE_FAN);
      glVertex3f(0., 0., -_height);
      for(short i = _segments; i >= 0; --i) {
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex2f(x, y);
      }
      glEnd();
    }


    ////////////////////////////////////////////////////////////////////////////
    /// Draw a wire cone with the base centered at the origin and tip pointed away
    /// from the camera.
    /// @param[in] _radius The radius of the base.
    /// @param[in] _height The height of the cone.
    /// @param[in] _segments The number of segments to use for the sides.
    inline void
    cone_frame(const GLfloat _radius, const GLfloat _height,
        const unsigned short _segments = 20)
    {
      circle_frame(_radius, _segments);

      GLfloat incr = 2 * PI / _segments;
      GLfloat x, y;

      glBegin(GL_LINES);
      for(short i = _segments; i > 0; --i) {
        glVertex3f(0., 0., -_height);
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex2f(x, y);
      }
      glEnd();
    }


    ////////////////////////////////////////////////////////////////////////////
    /// Draw a cylinder centered at the origin and oriented along the z-axis.
    /// @param[in] _radius The cylinder radius.
    /// @param[in] _length The length perpendicular to the radius.
    /// @param[in] _segments The number of segments to use for the side wall.
    inline void
    cylinder(const GLfloat _radius, const GLfloat _length,
        const unsigned short _segments = 20)
    {
      GLfloat halfLength = _length / 2.;

      // Draw front cap.
      glTranslatef(0., 0., halfLength);
      circle(_radius, _segments);

      // Draw back cap.
      glTranslatef(0., 0., -_length);
      glRotatef(180, 0, 1, 0);
      circle(_radius, _segments);
      glRotatef(-180, 0, 1, 0);

      // Draw walls.
      glTranslatef(0., 0., halfLength);
      GLfloat incr = 2 * PI / _segments;
      GLfloat x, y;

      glBegin(GL_TRIANGLE_STRIP);
      for(short i = 0; i <= _segments; ++i) {
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex3f(x, y,  halfLength);
        glVertex3f(x, y, -halfLength);
      }
      glEnd();
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Draw a wire cylinder.
    inline void
    cylinder_frame(const GLfloat _radius, const GLfloat _length,
        const unsigned short _segments = 20)
    {
      GLfloat halfLength = _length / 2.;

      // Draw front cap.
      glTranslatef(0., 0., halfLength);
      circle_frame(_radius, _segments);

      // Draw back cap.
      glTranslatef(0., 0., -_length);
      circle_frame(_radius, _segments);

      // Draw walls.
      glTranslatef(0., 0., halfLength);
      GLfloat incr = 2 * PI / _segments;
      GLfloat x, y;

      glBegin(GL_LINES);
      for(short i = 0; i < _segments; ++i) {
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex3f(x, y,  halfLength);
        glVertex3f(x, y, -halfLength);
      }
      glEnd();
    }

    ///@}

  }

}

#endif
