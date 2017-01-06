#include "gl_widget.h"

#include <array>

#include <QtGui>

#include "glutils/camera.h"
#include "glutils/color.h"

#include "base_visualization.h"
#include "std_input_manager.h"


using namespace std;
using namespace glutils;

/*------------------------------- Construction -------------------------------*/

gl_widget::
gl_widget(QWidget* _parent) :
    QGLWidget(_parent)
{
  // Set widget size and initialize GL.
  setMinimumSize(800, 600);
  makeCurrent();
  initializeGL();

  // Create camera and set default position.
  m_camera = new glutils::camera();
  m_camera->position(glutils::vector3f{0,0,10}, glutils::vector3f{0,0,0},
      glutils::vector3f{0,1,0});

  // Create controller manager.
  m_input_manager = new std_input_manager(this);

  // Set up rendering clock.
  m_clock = new QTimer(this);
  m_clock->setInterval(33);
  connect(m_clock, SIGNAL(timeout()),
          this,      SLOT(updateGL()));

  // Use the mouse and refrain from propagating unhandled events to (_parent).
  this->setMouseTracking(true);
  this->setAttribute(Qt::WA_NoMousePropagation);

  // Set focus policy to accept key press events.
  this->setFocusPolicy(Qt::StrongFocus);
}


gl_widget::
~gl_widget()
{
  delete m_camera;
}

/*------------------------------ Accessors -----------------------------------*/

glutils::camera*
gl_widget::
camera() const
{
  return m_camera;
}


base_visualization*
gl_widget::
visualization() const
{
  return m_visualization;
}


void
gl_widget::
visualization(base_visualization* const _s)
{
  m_visualization = _s;
}

/*--------------------------- Simulation Control -----------------------------*/

void
gl_widget::
start()
{
  emit status_message(QString("Starting visualization..."), 2000);

  m_visualization->start();
  m_clock->start();
}


void
gl_widget::
reset()
{
  emit status_message(QString("Stopping visualization..."), 2000);

  m_clock->stop();
  m_visualization->reset();
  m_camera->position(glutils::vector3f{0,0,10}, glutils::vector3f{0,0,0},
      glutils::vector3f{0,1,0});
  updateGL();
}


/*-------------------------------- Selection ---------------------------------*/

/// Make a box from two QPoints. Returns x-center, y-center, width, height.
inline
std::array<size_t, 4>
make_box(const QPoint& _p1, const QPoint& _p2)
{
  const size_t x = static_cast<size_t>((_p1.x() + _p2.x()) / 2);
  const size_t y = static_cast<size_t>((_p1.y() + _p2.y()) / 2);
  const size_t w = static_cast<size_t>(std::abs(_p1.x() - _p2.x()));
  const size_t h = static_cast<size_t>(std::abs(_p1.y() - _p2.y()));
  return std::array<size_t, 4>{x, y, w, h};
}


void
gl_widget::
select(const size_t _x, const size_t _y, const size_t _w, const size_t _h)
{
  // Qt convention is backward from OpenGL. Invert y coordinate to fix.
  m_visualization->render_select(_x, QWidget::height() - _y, _w, _h);
}


void
gl_widget::
select(QPoint _p1, QPoint _p2)
{
  auto b = make_box(_p1, _p2);
  select(b[0], b[1], b[2], b[3]);
}


void
gl_widget::
select(QPoint _p)
{
  select(_p.rx(), _p.ry(), 1, 1);
}


void
gl_widget::
hover(const size_t _x, const size_t _y, const size_t _w, const size_t _h)
{
  // Qt convention is backward from OpenGL. Invert y coordinate to fix.
  m_visualization->render_hover(_x, QWidget::height() - _y, _w, _h);
}


void
gl_widget::
hover(QPoint _p1, QPoint _p2)
{
  auto b = make_box(_p1, _p2);
  hover(b[0], b[1], b[2], b[3]);
}


void
gl_widget::
hover(QPoint _p)
{
  hover(_p.rx(), _p.ry(), 1, 1);
}

/*------------------------------- GL Functions -------------------------------*/

void
gl_widget::
initializeGL()
{
  // Set clear color to black
  glClearColor(0., 0., 0., 1.);

  // Enable material coloring with glColor
  glEnable(GL_COLOR_MATERIAL);

  // Create some lights to show 3d objects as 3d
  glEnable(GL_LIGHTING);
  glLightfv(GL_LIGHT0, GL_AMBIENT, color::black);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, color::white);
  glLightfv(GL_LIGHT0, GL_SPECULAR, color::white);
  glLightfv(GL_LIGHT1, GL_AMBIENT, color::black);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, color::white);
  glLightfv(GL_LIGHT1, GL_SPECULAR, color::white);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);

  // Position the lights
  GLfloat light0Pos[] = {250., 250., 250., 1.};
  GLfloat light1Pos[] = {-250., 250., -250., 1.};
  glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
  glLightfv(GL_LIGHT1, GL_POSITION, light1Pos);

  // Set material response to light
  glMaterialfv(GL_FRONT, GL_AMBIENT, color::light_grey);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, color::light_grey);
  glMaterialfv(GL_FRONT, GL_SPECULAR, color::white);
  glMaterialf(GL_FRONT, GL_SHININESS, 100.);

  // Enable depth test
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glClearDepth(1.0);

  // Enable alpha blending
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Enable point anti-aliasing
  glEnable(GL_POINT_SMOOTH);

  // Enable smooth shading
  glShadeModel(GL_SMOOTH);

  // Enable back-face culling to see through boundaries
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
}


void
gl_widget::
paintGL()
{
  // Erase the previous frame.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set view and lighting.
  glLoadIdentity();
  m_camera->apply_view();

  // Render the scene.
  if(m_visualization)
    m_visualization->render();

  // Push data.
  glFlush();
}


void
gl_widget::
resizeGL(int _w, int _h)
{
  // Set the GL window size
  glViewport(0, 0, GLint(_w), GLint(_h));

  // Set the viewing volume
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-_w / 600., _w / 600., -_h / 600., _h / 600., 1, 10000);

  // Return to model view mode
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

/*------------------------------- Input Forwarding ---------------------------*/

void
gl_widget::
mousePressEvent(QMouseEvent* _e)
{
  m_input_manager->mousePressEvent(_e);
}


void
gl_widget::
mouseReleaseEvent(QMouseEvent* _e)
{
  m_input_manager->mouseReleaseEvent(_e);
}


void
gl_widget::
mouseMoveEvent(QMouseEvent* _e)
{
  m_input_manager->mouseMoveEvent(_e);
}


void
gl_widget::
mouseDoubleClickEvent(QMouseEvent* _e)
{
  m_input_manager->mouseDoubleClickEvent(_e);
}


void
gl_widget::
wheelEvent(QWheelEvent* _e)
{
  m_input_manager->wheelEvent(_e);
}


void
gl_widget::
keyPressEvent(QKeyEvent* _e)
{
  m_input_manager->keyPressEvent(_e);
}


void
gl_widget::
keyReleaseEvent(QKeyEvent* _e)
{
  m_input_manager->keyReleaseEvent(_e);
}


void
gl_widget::
tabletEvent(QTabletEvent* _e)
{
  m_input_manager->tabletEvent(_e);
}

/*----------------------------------------------------------------------------*/
