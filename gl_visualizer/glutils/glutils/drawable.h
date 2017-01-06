#ifndef GLUTILS_DRAWABLE_H_
#define GLUTILS_DRAWABLE_H_

#include <atomic>
#include <list>
#include <mutex>
#include <queue>

#include "glutils/color.h"
#include "glutils/gltraits.h"


namespace glutils {

  //////////////////////////////////////////////////////////////////////////////
  /// An abstract interface class for drawable objects.
  ///
  /// @details Derived classes must implement the build function to describe the
  ///          commands needed to draw the object. Selection and highlighting
  ///          lists are optional. The object maintains a queue of transforms
  ///          that allow one to store multiple pre-computed positions. On each
  ///          render call, the front of the transform queue will be used for
  ///          the modelview matrix. Use update_transform to move on to the next
  ///          queued position (if available).
  //////////////////////////////////////////////////////////////////////////////
  class drawable
  {

    private:

      ///@name Internal State
      ///@{

      const GLuint m_selection_id;            ///< ID for selection rendering.
      const color m_picking_color;            ///< Color for color-picking.

      std::atomic<bool> m_selected{false};    ///< Is selected?
      std::atomic<bool> m_highlighted{false}; ///< Is highlighted?

      /// The transform queue.
      std::queue<transform, std::list<transform>> m_transforms;

    protected:

      bool m_initialized{false};              ///< Is initialized?

      ///@}

    public:

      ///@name Construction
      ///@{

      drawable() noexcept;
      drawable(const drawable& _d) noexcept;
      drawable(drawable&& _d) noexcept;

      virtual ~drawable() = default;

      ///@}
      ///@name Rendering
      ///@{

      /// Render this object in the scene.
      void render();

      /// Render this object for gl selection.
      void render_select();

      /// Render this object for color picking.
      void render_color_pick();

      ///@}
      ///@name Transform
      ///@{
      /// Drawables store a queue of transformations, which allows
      /// precomputation of future rendering positions.

      /// Push a new transform onto the queue.
      void push_transform(const transform& _t) noexcept;

      /// Update to the next queued transform. If no queued transform is
      /// available, the current one will be retained.
      void update_transform() noexcept;

      /// Replace the current transform queue with a single identity transform.
      void clear_transform() noexcept;

      ///@}
      ///@name Highlighting and Selection
      ///@{
      /// Two sets of decorations are supported for selected and highlighted
      /// objects. These functions control whether those decorations are rendered.
      /// These functions can be overriden to allow base classes to react to
      /// changes in highlight/selection. If that is done, the corresponding
      /// base class method here should also be called to ensure the rendering
      /// effect occurs.

      /// Get this object's selection name. Should only be needed by selector.
      GLuint id() const noexcept;

      /// Get this object's picking color. Should only be needed by color
      /// picker.
      const color& picking_color() const noexcept;


      /// Is this object 'selected'?
      bool selected() const noexcept;

      /// Enable drawing this object's 'selected' decorations.
      virtual void select() noexcept;

      /// Disable drawing this object's 'selected' decorations.
      virtual void deselect() noexcept;


      /// Is this object 'highlighted'?
      bool highlighted() const noexcept;

      /// Enable drawing this object's 'highlighted' components.
      virtual void highlight() noexcept;

      /// Disable drawing this object's 'highlighted' components.
      virtual void unhighlight() noexcept;

      ///@}

    protected:

      ///@name Drawing Instructions
      ///@{
      /// These functions specify how to draw this object. In addition to the
      /// base graphic, decorations for selected and highlighted objects are
      /// (optionally) supported. The modelview matrix will be saved/restored
      /// prior to/after drawing, so there is no need to do so explicitly.

      virtual void initialize(); ///< Init to perform on first draw.

      virtual void draw() = 0;           ///< Render the object.
      virtual void draw_select() = 0;    ///< Render for selection (no colors).
      virtual void draw_selected() {}    ///< Render the selection decorations.
      virtual void draw_highlighted() {} ///< Render the highlight decorations.

      ///@}

    private:

      ///@name Selection ID Generation
      ///@{

      static std::mutex m_selection_id_gate; ///< Lock for the id function.

      /// Generate a unique selection id.
      static GLuint generate_selection_id() noexcept;

      /// Generate a unique picking color.
      static color generate_picking_color() noexcept;

      ///@}

  };

}

#endif
