#ifndef GLUTILS_DRAWABLE_CALL_LIST_H_
#define GLUTILS_DRAWABLE_CALL_LIST_H_

#include <memory>

#include "glutils/drawable.h"
#include "glutils/gltraits.h"

namespace glutils {

  class call_list_set;

  //////////////////////////////////////////////////////////////////////////////
  /// A drawable object that can be represented with a call list.
  //////////////////////////////////////////////////////////////////////////////
  class drawable_call_list :
      public drawable
  {

    ///@name Local Types
    ///@{

    typedef std::shared_ptr<call_list_set> call_list_ptr;

    ///@}
    ///@name Internal State
    ///@{

    call_list_ptr m_lists; ///< Shared pointer to a set of call lists.

    ///@}

    public:

      ///@name Construction
      ///@{

      /// Create a new drawable with a new call list.
      drawable_call_list();

      /// Create a drawable from an existing call list.
      /// @param[in] _cl A shared pointer to an existing call list set, which
      ///                will be used in place of generating a new call list.
      drawable_call_list(call_list_ptr _cl) noexcept;

      virtual void rebuild() override;

      virtual ~drawable_call_list();

      ///@}

    protected:

      ///@name Call List Specification
      ///@{
      /// These functions should be defined by subclasses to specify the
      /// appropriate rendering instructions.

      /// Instructions for drawing this object at the origin, from the standard
      /// OpenGL perspective.
      virtual void build() = 0;

      /// As build, but purely geometry only (no colors, textures, etc).
      virtual void build_select() = 0;

      /// Additional instructions for drawing the selection decorations on the
      /// object. These will be drawn over the object when it is selected.
      virtual void build_selected() {}

      /// Additional instructions for drawing the highlight decorations on the
      /// object. These will be drawn over the object when it is highlighted,
      /// after any decorations for selection.
      virtual void build_highlighted() {}

    private:

      ///@}
      ///@name drawable Overrides
      ///@{

      virtual void initialize() override final; ///< Compile the call lists.

      virtual void draw() override final;
      virtual void draw_select() override final;
      virtual void draw_selected() override final;
      virtual void draw_highlighted() override final;

      ///@}
  };

}

#endif
