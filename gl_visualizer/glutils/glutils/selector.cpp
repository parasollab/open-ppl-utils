#include "glutils/selector.h"

#include "glutils/drawable.h"

#include <algorithm>
#include <utility>

#include <GL/glu.h>

#include <iostream>

namespace glutils {

  /*-------------------------- Map Management --------------------------------*/

  void
  selector::
  add_drawable(drawable* _d)
  {
    m_drawables[_d->id()] = _d;
  }


  void
  selector::
  remove_drawable(drawable* _d)
  {
    auto iter = m_drawables.find(_d->id());
    if(iter != m_drawables.end()) {
      m_hits.erase(iter->first);
      m_unhits.erase(iter->first);
      m_drawables.erase(iter);
    }
  }

  /*------------------------------ Selection ---------------------------------*/

  std::vector<drawable*>
  selector::
  hits() const noexcept
  {
    std::vector<drawable*> hits;
    hits.reserve(m_hits.size());
    for(const auto& hit : m_hits)
      hits.push_back(m_drawables.at(hit));
    return hits;
  }


  std::vector<drawable*>
  selector::
  unhits() const noexcept
  {
    std::vector<drawable*> unhits;
    unhits.reserve(m_unhits.size());
    for(const auto& unhit : m_unhits)
      unhits.push_back(m_drawables.at(unhit));
    return unhits;
  }


  void
  selector::
  select(size_t _x, size_t _y, size_t _w, size_t _h)
  {
    if(m_debug)
      std::cout << "glutils::selector: "
                << "Selecting at (" << _x << ", " << _y << ") with range ["
                << _w << ":" << _h << "]" << std::endl;

    // Set the selection buffer, switch to selection mode, and initalize the
    // name stack.
    glSelectBuffer(m_buffer_size, m_buffer);
    glRenderMode(GL_SELECT);
    glInitNames();

    // Resize the projection matrix (i.e. picking window) to select only items
    // near the target.
    // Grab current viewport.
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    // Grab current projection matrix.
    GLdouble projection[16];
    glGetDoublev(GL_PROJECTION_MATRIX, projection);

    // Save the current projection matrix.
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();

    glLoadIdentity();                        // Start from origin.
    gluPickMatrix(_x, _y, _w, _h, viewport); // Set picking window.
    glMultMatrixd(projection);               // Set current view.

    // Render each drawable for selection.
    glMatrixMode(GL_MODELVIEW);
    for(auto& d : m_drawables)
      d.second->render_select();

    // Restore the previous projection matrix and return to modelview mode.
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);

    // Select everything in the target window if it was larger than 5 pixels.
    const bool select_all = _w > 5 or _h > 5;

    // Return to render mode and parse the hit buffer.
    const GLint num_hits = glRenderMode(GL_RENDER);
    parse_hit_buffer(num_hits, select_all);
  }

  /*--------------------------------- Debug ----------------------------------*/

  void
  selector::
  debug(const bool _show)
  {
    m_debug = _show;
  }

  /*-------------------------------- Helpers ---------------------------------*/

  void
  selector::
  parse_hit_buffer(const GLint _num_hits, const bool _all)
  {
    if(m_debug)
      std::cout << "glutils::selector: "
                << "Parsing " << (_all ? "all" : "nearest") << " of "
                << _num_hits << " hits." << std::endl;

    std::pair<float, GLuint> nearest{std::numeric_limits<float>::max(), 0};
    hit_list hits;

    // Parse current hits.
    GLuint* buffer = m_buffer - 1;
    for(GLint i = 0; i < _num_hits; ++i) {
      // Get the number of names on the stack for this hit.
      GLuint num_names = *++buffer;

      // Get the near distance for this hit.
      float z_near = static_cast<float>(*++buffer) / 0x7fffffff;
      // float z_far  = static_cast<float>(m_buffer[2]) / 0x7fffffff;
      ++buffer; // skip far distance.

      // Add each name to the hit list.
      for(size_t n = 0; n < num_names; ++n) {
        // Add this name.
        const GLuint name = *++buffer;
        hits.insert(name);

        // Track nearest hit.
        if(z_near < nearest.first) {
          nearest.first = z_near;
          nearest.second = name;
        }
      }
    }

    // If we just want the nearest, mark every other previous hit as unhit.
    if(!_all) {
      m_unhits = std::move(m_hits);
      m_unhits.erase(nearest.second);
      m_hits.clear();
      if(!hits.empty())
        m_hits.insert(nearest.second);
    }
    // Otherwise, sort out which names are hits and unhits.
    else {
      // The new hits are those that are now selected, but previously weren't.
      // The new unhits are that were previously selected, but aren't anymore.
      hit_list new_hits, new_unhits;

      std::set_difference(hits.begin(), hits.end(), m_hits.begin(), m_hits.end(),
          std::inserter(new_hits, new_hits.begin()));
      std::set_difference(m_hits.begin(), m_hits.end(), hits.begin(), hits.end(),
          std::inserter(new_unhits, new_unhits.begin()));

      m_unhits = std::move(new_unhits);
      m_hits = std::move(new_hits);
    }

    if(m_debug) {
      std::cout << "\tHits:";
      for(const auto& hit : m_hits)
        std::cout << " " << hit;

      std::cout << "\n\tUnhits:";
      for(const auto& unhit : m_unhits)
        std::cout << " " << unhit;

      std::cout << std::endl;
    }
  }

  /*--------------------------------------------------------------------------*/

}
