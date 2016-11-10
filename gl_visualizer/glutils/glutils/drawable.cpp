#include "glutils/drawable.h"

namespace glutils {

  /*------------------------------ Construction ------------------------------*/

  drawable::
  drawable() noexcept :
      m_selection_id(generate_selection_id())
  {
    push_transform(identity_transform());
  }


  drawable::
  drawable(const drawable& _d) noexcept :
      m_selection_id(generate_selection_id())
  {
    push_transform(identity_transform());
  }


  drawable::
  drawable(drawable&& _d) noexcept :
      m_selection_id(_d.m_selection_id),
      m_selected(_d.m_selected.load()),
      m_highlighted(_d.m_highlighted.load()),
      m_transforms(std::move(_d.m_transforms))
  {
    _d.m_selected = false;
    _d.m_highlighted = false;
  }

  /*-------------------------------- Rendering -------------------------------*/

  void
  drawable::
  render()
  {
    glPushMatrix();

    if(!m_initialized)
      initialize();

    apply_transform(m_transforms.front());

    draw();
    if(m_selected)
      draw_selected();
    if(m_highlighted)
      draw_highlighted();;

    glPopMatrix();
  }


  void
  drawable::
  render_select()
  {
    glPushMatrix();

    apply_transform(m_transforms.front());

    glPushName(m_selection_id);
    draw();
    glPopName();

    glPopMatrix();
  }

  /*------------------------------- Transform --------------------------------*/

  void
  drawable::
  push_transform(const transform& _t) noexcept
  {
    m_transforms.push(_t);
  }


  void
  drawable::
  update_transform() noexcept
  {
    if(m_transforms.size() > 1)
      m_transforms.pop();
  }

  /*-------------------------- Drawing Instructions --------------------------*/

  void
  drawable::
  initialize()
  {
    m_initialized = true;
  }

  /*------------------------ Selection ID Generation -------------------------*/

  std::mutex
  drawable::m_selection_id_gate;


  const GLuint
  drawable::
  generate_selection_id() noexcept
  {
    m_selection_id_gate.lock();
    static GLuint next = 0;
    GLuint id = ++next;
    m_selection_id_gate.unlock();
    return id;
  }

  /*--------------------------------------------------------------------------*/

}
