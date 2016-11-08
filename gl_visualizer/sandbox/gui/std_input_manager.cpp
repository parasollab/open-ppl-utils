#include "std_input_manager.h"

#include "gl_widget.h"

#include "nonstd/numerics.h"

#include "glutils/camera.h"
#include "glutils/draw.h"


/*------------------------------ Button State --------------------------------*/

void
std_input_manager::
button_state::
reset() noexcept
{
  drag = 0;
  click.setX(0);
  click.setY(0);
}

/*------------------------------ Construction --------------------------------*/

std_input_manager::
std_input_manager(gl_widget* _gl) :
    QWidget(_gl), m_gl(_gl)
{
  connect(this, SIGNAL(mouse_selection(QPoint)),
          m_gl, SLOT(select(QPoint)));
  connect(this, SIGNAL(mouse_selection(QPoint, QPoint)),
          m_gl, SLOT(select(QPoint, QPoint)));

  connect(this, SIGNAL(mouse_hover(QPoint)),
          m_gl, SLOT(hover(QPoint)));
  connect(this, SIGNAL(mouse_hover(QPoint, QPoint)),
          m_gl, SLOT(hover(QPoint, QPoint)));
}

/*------------------------------ Input Events --------------------------------*/

void
std_input_manager::
mousePressEvent(QMouseEvent* _e)
{
  switch(_e->button()) {
    case Qt::LeftButton:
      click_left(_e, true);
      break;
    case Qt::MiddleButton:
      click_middle(_e, true);
      break;
    case Qt::RightButton:
      click_right(_e, true);
      break;
    default:;
  }
}


void
std_input_manager::
mouseReleaseEvent(QMouseEvent* _e)
{
  switch(_e->button()) {
    case Qt::LeftButton:
      click_left(_e, false);
      break;
    case Qt::MiddleButton:
      click_middle(_e, false);
      break;
    case Qt::RightButton:
      click_right(_e, false);
      break;
    default:;
  }
}


void
std_input_manager::
mouseMoveEvent(QMouseEvent* _e)
{
  if(!_e->buttons()) {
    // If no buttons were pressed during the movement, just update mouse position.
    emit mouse_hover(_e->pos());
  }
  else {
    // There is at least one button held down. Determine the displacement and
    // whether the deadzone was exceeded.
    const QPoint delta = _e->pos() - m_hover;
    const bool deadzone_exceeded = delta.manhattanLength() >= m_deadzone;

    // Process the drag even for each button.
    if(_e->buttons() & Qt::LeftButton && (m_left.drag || deadzone_exceeded))
      drag_left(_e, delta);
    if(_e->buttons() & Qt::MiddleButton && (m_middle.drag || deadzone_exceeded))
      drag_middle(_e, delta);
    if(_e->buttons() & Qt::RightButton && (m_right.drag || deadzone_exceeded))
      drag_right(_e, delta);
  }

  // Update last hover position.
  m_hover = _e->pos();
}


void
std_input_manager::
mouseDoubleClickEvent(QMouseEvent*)
{ }


void
std_input_manager::
wheelEvent(QWheelEvent* _e)
{
  // Scale wheel distance quadratically to make vigorous scrolling zoom faster.
  int wheelDist = _e->delta() / 60.;
  wheelDist *= wheelDist * nonstd::sign(wheelDist);
  m_gl->camera()->zoom(wheelDist);
}


void
std_input_manager::
keyPressEvent(QKeyEvent*)
{ }


void
std_input_manager::
keyReleaseEvent(QKeyEvent*)
{ }


void
std_input_manager::
tabletEvent(QTabletEvent*)
{ }

/*----------------------------------------------------------------------------*/

void
std_input_manager::
click_left(QMouseEvent* _e, const bool _press)
{
  if(_press) {
    // Handle click.
    m_left.click = _e->pos();
  }
  else {
    // Handle unclick.
    // If there is an active middle or right drag, do nothing.
    if(m_right.drag || m_middle.drag)
      ;
    // If there is no active left drag, this is a point selection.
    else if(!m_left.drag)
      emit mouse_selection(_e->pos());
    // Otherwise, this is an area selection.
    else
      emit mouse_selection(m_left.click, _e->pos());

    m_left.reset();
  }
}


void
std_input_manager::
click_middle(QMouseEvent* _e, const bool _press)
{
  if(_press) {
    // Handle click.
    m_middle.click = _e->pos();
  }
  else {
    // Handle unclick.
    m_middle.reset();
  }
}


void
std_input_manager::
click_right(QMouseEvent* _e, const bool _press)
{
  if(_press) {
    // Handle click.
    m_right.click = _e->pos();
  }
  else {
    // Handle unclick.
    m_right.reset();
  }
}


void
std_input_manager::
drag_left(QMouseEvent* _e, const QPoint& _delta)
{
  m_left.drag += _delta.manhattanLength();

  // If the right mouse button is also down, do nothing. Otherwise, emit hover
  // signal for the box from the click point to here.
  if(!(_e->buttons() & Qt::RightButton))
    emit mouse_hover(m_left.click, _e->pos());
}


void
std_input_manager::
drag_middle(QMouseEvent* _e, const QPoint& _delta)
{
  m_middle.drag += _delta.manhattanLength();

  // Pan the camera view.
  auto camera = m_gl->camera();
  camera->pan(_delta.x() * m_sensitivity, -_delta.y() * m_sensitivity);
}


void
std_input_manager::
drag_right(QMouseEvent* _e, const QPoint& _delta)
{
  m_right.drag += _delta.manhattanLength();

  // If the left mouse button is also held down, rotate the camera about itself.
  // Otherwise, orbit the camera about the origin.
  auto camera = m_gl->camera();
  if(_e->buttons() & Qt::LeftButton) {
    camera->rotate(-_delta.x() * m_sensitivity * glutils::RadPerDeg, camera->y());
    camera->rotate(-_delta.y() * m_sensitivity * glutils::RadPerDeg, camera->x());
  }
  else {
    camera->rotate(-_delta.x() * m_sensitivity * glutils::RadPerDeg,
        camera->y(), glutils::vector3f());
    camera->rotate(-_delta.y() * m_sensitivity * glutils::RadPerDeg,
        camera->x(), glutils::vector3f());
  }
}

/*----------------------------------------------------------------------------*/
