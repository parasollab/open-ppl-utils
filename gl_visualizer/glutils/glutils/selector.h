#ifndef GLUTILS_SELECTOR_H_
#define GLUTILS_SELECTOR_H_

#include "glutils/gltraits.h"

#include <cstddef>
#include <set>
#include <unordered_map>
#include <vector>

namespace glutils {

  class drawable;

  //////////////////////////////////////////////////////////////////////////////
  /// A gl selection mapper for instances of drawable.
  //////////////////////////////////////////////////////////////////////////////
  class selector final
  {

    ///@name Local Types
    ///@{

    /// A mapping from selection id to drawable.
    typedef std::unordered_map<GLuint, drawable*> selection_map;

    /// A list of selected ids.
    typedef std::set<GLuint> hit_list;

    ///@}
    ///@name Internal State
    ///@{

    selection_map m_drawables; ///< The selection map.
    hit_list m_hits;           ///< The id's that were selected on last check.
    hit_list m_unhits;         ///< The id's that were deselected by last check.

    static constexpr size_t m_buffer_size{1024}; ///< The size of the hit buffer.
    GLuint m_buffer[m_buffer_size];              ///< This selector's hit buffer.

    bool m_debug{false};       ///< Show debugging messages?

    ///@}

    public:

      ///@}
      ///@name Map Management
      ///@{

      //////////////////////////////////////////////////////////////////////////
      /// Add a drawable to the selection map.
      /// @param _d The drawable to add.
      void add_drawable(drawable* _d);

      //////////////////////////////////////////////////////////////////////////
      /// Remove a drawable from the selection map.
      /// @param _d The drawable to remove.
      void remove_drawable(drawable* _d);

      ///@}
      ///@name Selection
      ///@{

      //////////////////////////////////////////////////////////////////////////
      /// Get the models that were selected on last attempt.
      std::vector<drawable*> hits() const noexcept;

      //////////////////////////////////////////////////////////////////////////
      /// Get the models that were unselected on last attempt.
      std::vector<drawable*> unhits() const noexcept;

      //////////////////////////////////////////////////////////////////////////
      /// Perform GL selection, marking the selected drawables as such. This
      /// must be called from the rendering thread.
      /// @param _x The X window coordinate for the pick box center. This uses
      ///           the GL convention, where X = 0 is the left of the screen.
      /// @param _y The Y window coordinate for the pick box center. This uses
      ///           the GL convention, where Y = 0 is the bottom of the screen.
      /// @param _w The pick box width.
      /// @param _h The pick box height.
      void select(size_t _x, size_t _y, size_t _w = 5, size_t _h = 5);

      ///@}
      ///@name Debug
      ///@{

      void debug(const bool _show); ///< Enable or disable debugging messages.

      ///@}

    private:

      ///@name Helpers
      ///@{

      //////////////////////////////////////////////////////////////////////////
      /// Parse the hit buffer and mark selected items.
      /// @param _num_hits The number of items in the hit buffer.
      /// @param _all Select all hits or just the closest?
      void parse_hit_buffer(const GLint _num_hits, const bool _all);

      ///@}

  };

}

#endif
