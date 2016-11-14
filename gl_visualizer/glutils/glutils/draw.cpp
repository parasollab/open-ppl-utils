#include "glutils/draw.h"


namespace glutils {

  namespace draw {

    /*--------------------------- 2D Primitives ------------------------------*/

    void
    circle(const GLfloat _radius, const size_t _segments)
    {
      GLfloat incr = 2. * PI / _segments;
      GLfloat x, y;

      glBegin(GL_TRIANGLE_FAN);
      glVertex2f(0., 0.);
      for(size_t i = 0; i <= _segments; ++i) {
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex2f(x, y);
      }
      glEnd();
    }


    void
    circle_frame(const GLfloat _radius, const size_t _segments)
    {
      GLfloat incr = 2 * PI / _segments;
      GLfloat x, y;

      glBegin(GL_LINE_LOOP);
      for(size_t i = 0; i <= _segments; ++i) {
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex2f(x, y);
      }
      glEnd();
    }

    /*--------------------------- 3D Primitives ------------------------------*/

    void
    sphere(const GLfloat _radius, const size_t _segments)
    {
      GLfloat oIncr = 2 * PI / _segments; // Angle increment for x,y coords.
      GLfloat zIncr = PI / _segments;     // Angle increment for z coords.
      GLfloat x, y, z, r;

      // Draw +zHat cap.
      glBegin(GL_TRIANGLE_FAN);
      glVertex3f(0, 0, _radius); // The +zHat pole.
      z = _radius * cos(zIncr);
      r = _radius * sin(zIncr);
      for(size_t i = 0; i <= _segments; ++i) {
        x = r * cos(oIncr * i);
        y = r * sin(oIncr * i);
        glVertex3f(x, y, z);
      }
      glEnd();

      // Draw main surface.
      GLfloat z2, r2;
      for(size_t j = 1; j < _segments; ++j) {
        glBegin(GL_TRIANGLE_STRIP);
        z  = _radius * cos(zIncr * j);
        r  = _radius * sin(zIncr * j);
        z2 = _radius * cos(zIncr * (j + 1));
        r2 = _radius * sin(zIncr * (j + 1));
        for(size_t i = 0; i <= _segments; ++i) {
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
      for(int i = static_cast<int>(_segments); i >= 0; --i) {
        x = r * cos(oIncr * i);
        y = r * sin(oIncr * i);
        glVertex3f(x, y, z);
      }
      glEnd();
    }


    void
    sphere_frame(const GLfloat _radius, const size_t _segments)
    {
      GLfloat oIncr = 2 * PI / _segments; // Angle increment for x,y coords.
      GLfloat zIncr = PI / _segments;     // Angle increment for z coords.
      GLfloat x, y, z, r;

      // Draw latitude lines.
      glBegin(GL_LINES);

      // Draw +zHat cap.
      z = _radius * cos(zIncr);
      r = _radius * sin(zIncr);
      for(size_t i = 0; i < _segments; ++i) {
        glVertex3f(0, 0, _radius);
        x = r * cos(oIncr * i);
        y = r * sin(oIncr * i);
        glVertex3f(x, y, z);
      }

      // Draw main surface.
      GLfloat z2, r2;
      for(size_t j = 1; j < _segments; ++j) {
        z  = _radius * cos(zIncr * j);
        r  = _radius * sin(zIncr * j);
        z2 = _radius * cos(zIncr * (j + 1));
        r2 = _radius * sin(zIncr * (j + 1));
        for(size_t i = 0; i <= _segments; ++i) {
          x = cos(oIncr * i);
          y = sin(oIncr * i);
          glVertex3f(x * r , y * r , z);
          glVertex3f(x * r2, y * r2, z2);
        }
      }

      // Draw -zHat cap.
      z = _radius * cos(zIncr * (_segments - 1));
      r = _radius * sin(zIncr * (_segments - 1));
      for(size_t i = _segments; i > 0; --i) {
        x = r * cos(oIncr * i);
        y = r * sin(oIncr * i);
        glVertex3f(x, y, z);
        glVertex3f(0, 0, -_radius);
      }
      glEnd();

      // Draw longitude lines.
      for(size_t i = 1; i < _segments; ++i) {
        glPushMatrix();
        z = _radius * cos(zIncr * i);
        r = _radius * sin(zIncr * i);

        glTranslatef(0, 0, z);
        circle_frame(r, _segments);
        glPopMatrix();
      }
    }


    void
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


    void
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


    void
    cone(const GLfloat _radius, const GLfloat _height,
        const size_t _segments)
    {
      circle(_radius, _segments);

      GLfloat incr = 2 * PI / _segments;
      GLfloat x, y;

      glBegin(GL_TRIANGLE_FAN);
      glVertex3f(0., 0., -_height);
      for(int i = static_cast<int>(_segments); i >= 0; --i) {
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex2f(x, y);
      }
      glEnd();
    }


    void
    cone_frame(const GLfloat _radius, const GLfloat _height,
        const size_t _segments)
    {
      circle_frame(_radius, _segments);

      GLfloat incr = 2 * PI / _segments;
      GLfloat x, y;

      glBegin(GL_LINES);
      for(size_t i = _segments; i > 0; --i) {
        glVertex3f(0., 0., -_height);
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex2f(x, y);
      }
      glEnd();
    }


    void
    cylinder(const GLfloat _radius, const GLfloat _length,
        const size_t _segments)
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
      for(size_t i = 0; i <= _segments; ++i) {
        x = _radius * cos(incr * i);
        y = _radius * sin(incr * i);
        glVertex3f(x, y,  halfLength);
        glVertex3f(x, y, -halfLength);
      }
      glEnd();
    }


    inline void
    cylinder_frame(const GLfloat _radius, const GLfloat _length,
        const size_t _segments)
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
      for(size_t i = 0; i < _segments; ++i) {
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
