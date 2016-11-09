#ifndef BASE_VISUALIZATION_H_
#define BASE_VISUALIZATION_H_

#include <cstddef>
#include <vector>

#include "glutils/drawable.h" // Includes rather than forward-declarations to
#include "glutils/selector.h" // bring these components to any base class.


////////////////////////////////////////////////////////////////////////////////
/// The base interface for a visualization that runs in the sandbox. It has
/// default storage for drawables, and a default selector/highlighter. The
/// rendering functions can be overriden to add additional instructions, but the
/// base class methods should always be called as well to retain the default
/// functionality.
////////////////////////////////////////////////////////////////////////////////
struct base_visualization
{

  ///@name Internal State
  ///@{

  std::vector<glutils::drawable*> m_drawables; ///< The drawable objects.

  glutils::selector* m_selector{nullptr};      ///< Selection helper.
  glutils::selector* m_highlighter{nullptr};   ///< Highlight helper.

  ///@}
  ///@name Construction
  ///@{

  base_visualization();
  virtual ~base_visualization();

  ///@}
  ///@name Required Interface
  ///@{

  /// Update the transforms for all drawable objects. This is not done in
  /// render to give more flexibility for multi-threaded visualizations.
  virtual void update();

  /// Define any instructions needed to render this visualization.
  virtual void render();

  /// Define any instructions to be executed when a selection takes place.
  /// @param _x The x coordinate of the picking box's center.
  /// @param _y The y coordinate of the picking box's center.
  /// @param _w The picking box width.
  /// @param _h The picking box height.
  virtual void render_select(const size_t _x, const size_t _y, const size_t _w,
      const size_t _h);

  /// Define any instructions to be executed when a hover event takes place.
  /// @param _x The x coordinate of the picking box's center.
  /// @param _y The y coordinate of the picking box's center.
  /// @param _w The picking box width.
  /// @param _h The picking box height.
  virtual void render_hover(const size_t _x, const size_t _y, const size_t _w,
      const size_t _h);

  /// Define any instructions to be executed at the start of the visualization.
  virtual void start() {}

  /// Define how to reset the visualization to it's initial state.
  virtual void reset() {}

  ///@}

};

#endif
