#include "glutils/call_list_set.h"

namespace glutils {

  /*------------------------------ Construction ------------------------------*/

  call_list_set::
  call_list_set(const size_t _num) :
      m_num(_num), m_lists(_num, 0)
  { }


  call_list_set::
  ~call_list_set()
  {
    glDeleteLists(m_lists.front(), m_num);
  }

  /*------------------------------- Accessors --------------------------------*/

  void
  call_list_set::
  initialize()
  {
    // Generate a call list for rendering, selection, and highlighting.
    m_lists.front() = glGenLists(m_num);

    // If we got 0, there was a problem getting the lists.
    if(m_lists.front() == 0)
      throw std::runtime_error("glutils::call_list_set::initialize() error: "
          "could not allocate OpenGL call lists for the object.");

    // Pre-compute the call lists since we will call them all the time.
    for(GLuint i = 1; i < m_num; ++i)
      m_lists[i] = m_lists.front() + i;

    m_ready = true;
  }

  /*--------------------------------------------------------------------------*/

}
