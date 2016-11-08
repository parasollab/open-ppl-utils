#include "base_visualization.h"

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
render()
{
  for(const auto& d : m_drawables)
    d->render();
}


void
base_visualization::
render_select(const size_t _x, const size_t _y, const size_t _w, const size_t _h)
{
  m_selector->select(_x, _y, _w, _h);
  for(auto hit : m_selector->hits())
    hit->select();
  for(auto unhit : m_selector->unhits())
    unhit->deselect();
}


void
base_visualization::
render_hover(const size_t _x, const size_t _y, const size_t _w, const size_t _h)
{
  /// @TODO Hover rendering is presently disabled because the performance sucks.
  ///       Need to find a better way to achieve this.
  //m_highlighter->select(_x, _y, _w, _h);
  //for(auto hit : m_highlighter->hits())
  //  hit->highlight();
  //for(auto unhit : m_highlighter->unhits())
  //  unhit->unhighlight();
}

/*----------------------------------------------------------------------------*/
