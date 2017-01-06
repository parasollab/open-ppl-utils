#include "base_visualization.h"

#include <algorithm>

/*----------------------------- Construction ---------------------------------*/

base_visualization::
base_visualization()
  : m_selector(new glutils::selector),
    m_highlighter(new glutils::selector)
{ }


base_visualization::
~base_visualization()
{
  delete m_selector;
  delete m_highlighter;
  for(auto d : m_drawables)
    delete d;
}

/*------------------------- Visualization Interface --------------------------*/

void
base_visualization::
add_drawable(glutils::drawable* _d)
{
  m_drawables.push_back(_d);
  m_selector->add_drawable(_d);
  m_highlighter->add_drawable(_d);
}


void
base_visualization::
remove_drawable(glutils::drawable* _d)
{
  auto iter = std::find(m_drawables.begin(), m_drawables.end(), _d);
  if(iter == m_drawables.end())
    return;

  m_drawables.erase(iter);
  m_selector->remove_drawable(_d);
  m_highlighter->remove_drawable(_d);
  delete _d;
}


void
base_visualization::
update()
{
  for(auto& d : m_drawables)
    d->update_transform();
}


void
base_visualization::
render()
{
  for(const auto& d : m_drawables)
    d->render();
}


void
base_visualization::
render_select(const size_t _x, const size_t _y, const size_t _w, const size_t _h)
{
  // Use color-picking for point selection and gl selection for area selection.
  //if(_w < 2 && _h < 2)
  //  m_selector->color_pick(_x, _y);
  //else
  //  m_selector->select(_x, _y, _w, _h);
  m_selector->select(_x, _y, _w, _h);

  for(auto unhit : m_selector->unhits())
    unhit->deselect();
  for(auto hit : m_selector->hits())
    hit->select();
}


void
base_visualization::
//render_hover(const size_t, const size_t, const size_t, const size_t)
render_hover(const size_t _x, const size_t _y, const size_t, const size_t)
{
  /// @TODO Hover rendering is presently disabled because the performance sucks.
  ///       Need to find a better way to achieve this.
  m_highlighter->color_pick(_x, _y);
  for(auto unhit : m_highlighter->unhits())
    unhit->unhighlight();
  for(auto hit : m_highlighter->hits())
    hit->highlight();
}

/*----------------------------------------------------------------------------*/
