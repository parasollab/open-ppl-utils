#include "glutils/drawable_call_list.h"

#include "glutils/call_list_set.h"


namespace glutils {

  /*------------------------------ Construction ------------------------------*/

  drawable_call_list::
  drawable_call_list() :
      m_lists(new call_list_set(4))
  { }

  drawable_call_list::
  drawable_call_list(call_list_ptr _cl) noexcept :
      m_lists(_cl)
  { }

  /*-------------------------- drawable Overrides ----------------------------*/

  void
  drawable_call_list::
  draw()
  {
    glCallList((*m_lists)[0]);
  }


  void
  drawable_call_list::
  draw_select()
  {
    glCallList((*m_lists)[1]);
  }


  void
  drawable_call_list::
  draw_selected()
  {
    glCallList((*m_lists)[2]);
  }


  void
  drawable_call_list::
  draw_highlighted()
  {
    glCallList((*m_lists)[3]);
  }

  /*------------------------------- Helpers ----------------------------------*/

  void
  drawable_call_list::
  initialize()
  {
    drawable::initialize();

    // If the lists are already compiled, there is nothing to do.
    if(m_lists->ready())
      return;

    // Initialize the call lists.
    m_lists->initialize();

    // Ensure that we start from the identity matrix for each list and don't
    // disturb the stack in use.
    glPushMatrix();

    glLoadIdentity();
    glNewList((*m_lists)[0], GL_COMPILE);
    build();
    glEndList();

    glLoadIdentity();
    glNewList((*m_lists)[1], GL_COMPILE);
    build_select();
    glEndList();

    glLoadIdentity();
    glNewList((*m_lists)[2], GL_COMPILE);
    build_selected();
    glEndList();

    glLoadIdentity();
    glNewList((*m_lists)[3], GL_COMPILE);
    build_highlighted();
    glEndList();

    glPopMatrix();
  }

  /*--------------------------------------------------------------------------*/

}
