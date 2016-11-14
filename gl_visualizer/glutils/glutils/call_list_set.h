#ifndef GLUTILS_CALL_LIST_SET_H_
#define GLUTILS_CALL_LIST_SET_H_

#include <cstddef>
#include <vector>

#include "gltraits.h"

namespace glutils {

  //////////////////////////////////////////////////////////////////////////////
  /// Manages a set of consecutive GL call lists.
  ///
  /// This is primarily useful for sharing call-lists amongst several objects.
  /// Wrapping up one of these objects in a shared_ptr allows you to get a set
  /// of call lists that are shared and only deleted when the last client is
  /// destroyed.
  ///
  /// @warning These objects should be deleted by the rendering thread to be
  ///          sure that the GL system is ready to accept the delete request.
  //////////////////////////////////////////////////////////////////////////////
  class call_list_set final
  {

    ///@name Internal State
    ///@{

    const size_t m_num;          ///< The number of lists to maintain.
    std::vector<GLuint> m_lists; ///< The call list ids.
    bool m_ready{false};         ///< Are the call lists initialized?

    ///@}

    public:

      ///@name Construction
      ///@{

      //////////////////////////////////////////////////////////////////////////
      /// Create a set of call lists.
      /// @param _num The number of call lists to reserve.
      call_list_set(const size_t _num);

      //////////////////////////////////////////////////////////////////////////
      /// Call lists are released on destruction.
      ~call_list_set();

      ///@}
      ///@name Accecssors
      ///@{

      //////////////////////////////////////////////////////////////////////////
      /// Are the call lists initialized?
      bool ready() const noexcept;

      //////////////////////////////////////////////////////////////////////////
      /// Generate the call lists.
      ///
      /// @warning This isn't done during construction to allow for the
      ///          possibility of constructing the objects outside the rendering
      ///          thread.
      void initialize();

      //////////////////////////////////////////////////////////////////////////
      /// Get the number of stored call lists.
      size_t size() const noexcept;

      //////////////////////////////////////////////////////////////////////////
      /// Access a stored call list.
      ///
      /// @warning It is assumed that clients will NOT request a non-existant
      ///          index in order to benefit from noexcept.
      GLuint operator[](const size_t _i) const noexcept;

      ///@}

  };

  /*---------------------------- Inlined Accessors ---------------------------*/

  inline
  bool
  call_list_set::
  ready() const noexcept
  {
    return m_ready;
  }


  inline
  size_t
  call_list_set::
  size() const noexcept
  {
    return m_num;
  }


  inline
  GLuint
  call_list_set::
  operator[](const size_t _i) const noexcept
  {
    return m_lists[_i];
  }

  /*--------------------------------------------------------------------------*/

}

#endif
