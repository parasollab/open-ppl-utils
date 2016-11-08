#include "example_visualization.h"

#include "glutils/color.h"
#include "glutils/draw.h"
#include "glutils/drawable_call_list.h"


/// Example of a drawable object.
class drawable_sphere :
    public glutils::drawable_call_list
{

  virtual void build() override
  {
    glColor4f(.5, .1, .5, .5);
    glLineWidth(1);
    glutils::draw::sphere(5);
    glutils::draw::sphere_frame(5);
  }

  virtual void build_selected() override
  {
    glColor4fv(glutils::color::yellow);
    glLineWidth(4);
    glutils::draw::sphere_frame(5);
  }

  virtual void build_highlighted() override
  {
    glLineWidth(1);
    glColor4f(0, .1, .6, .4);
    glutils::draw::sphere(5.1);
  }

};

/*----------------------------- Construction ---------------------------------*/

example_visualization::
example_visualization()
{
  glutils::drawable* sphere = new drawable_sphere;
  m_drawables.push_back(sphere);

  sphere = new drawable_sphere;
  glutils::transform t;
  t.translation() = glutils::vector3f{-5, 0, 0};
  sphere->push_transform(t);
  sphere->update_transform();
  m_drawables.push_back(sphere);

  sphere = new drawable_sphere;
  t.translation() = glutils::vector3f{5, 0, 0};
  sphere->push_transform(t);
  sphere->update_transform();
  m_drawables.push_back(sphere);

  // Base class will delete the drawables on destruction.
  for(auto d : m_drawables) {
    m_selector->add_drawable(d);
    m_highlighter->add_drawable(d);
  }
}

/*----------------------------------------------------------------------------*/
