#ifndef GLUTILS_DRAWABLE_H_
#define GLUTILS_DRAWABLE_H_

#include <atomic>
#include <list>
#include <mutex>
#include <queue>

#include "glutils/gltraits.h"

namespace glutils {

  //////////////////////////////////////////////////////////////////////////////
  /// A base class for drawable objects. Derived classes must implement the
  /// build function to describe the commands needed to l call list for drawing
  /// the object. Selection and highlighting lists are optional.
  //////////////////////////////////////////////////////////////////////////////
  class drawable
  {

    private:

      ///@name Internal State
      ///@{

      const GLuint m_selection_id;            ///< ID for selection rendering.

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

      //////////////////////////////////////////////////////////////////////////
      /// Render this object in the scene.
      ///
      /// The default implementation is to render the call lists as appropriate.
      /// It can be overriden by derived classes to use other drawing methods.
      void render();

      //////////////////////////////////////////////////////////////////////////
      /// Render this object for gl selection.
      void render_select();

      ///@}
      ///@name Transform
      ///@{
      /// Drawables store a queue of transformations, which allows
      /// precomputation of future rendering positions.

      /// Push a new transform onto the queue.
      void push_transform(const transform& _t) noexcept;

      /// Update to the next queued transform.
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

      /// Is this object selected?
      bool selected() const noexcept;

      /// Enable drawing this object's 'selected' decorations.
      virtual void select() noexcept;

      /// Disable drawing this object's 'selected' decorations.
      virtual void deselect() noexcept;


      /// Is this object highlighted?
      bool highlighted() const noexcept;

      /// Enable drawing this object's 'highlighted' components.
      virtual void highlight() noexcept;

      /// Disable drawing this object's 'highlighted' components.
      virtual void unhighlight() noexcept;

      /// Get this object's selection name. Should only be needed by selector.
      GLuint id() const noexcept;

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
      virtual void draw_selected() {}    ///< Render the selection decorations.
      virtual void draw_highlighted() {} ///< Render the highlight decorations.

      ///@}

    private:

      ///@name Selection ID Generation
      ///@{

      static std::mutex m_selection_id_gate; ///< Lock for the id function.

      //////////////////////////////////////////////////////////////////////////
      /// Generate a unique selection id.
      static GLuint generate_selection_id() noexcept;

      ///@}

  };

  /*------------------- Inlined Highlighting and Selection -------------------*/

  inline
  bool
  drawable::
  selected() const noexcept
  {
    return m_selected.load();
  }


  inline
  void
  drawable::
  select() noexcept
  {
    m_selected = true;
  }


  inline
  void
  drawable::
  deselect() noexcept
  {
    m_selected = false;
  }


  inline
  bool
  drawable::
  highlighted() const noexcept
  {
    return m_highlighted.load();
  }


  inline
  void
  drawable::
  highlight() noexcept
  {
    m_highlighted = true;
  }


  inline
  void
  drawable::
  unhighlight() noexcept
  {
    m_highlighted = false;
  }


  inline
  GLuint
  drawable::
  id() const noexcept
  {
    return m_selection_id;
  }

  /*--------------------------------------------------------------------------*/

}

#endif
