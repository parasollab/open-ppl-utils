#ifndef GL_UTILS_TRIANGULATED_MODEL
#define GL_UTILS_TRIANGULATED_MODEL

#include <array>
#include <cstddef>
#include <vector>

#include "glutils/gltraits.h"

namespace glutils {

  //////////////////////////////////////////////////////////////////////////////
  /// A triangular, right-handed facet for a polygonal model.
  //////////////////////////////////////////////////////////////////////////////
  class triangle_facet {

    ///@name Local Types
    ///@{

    typedef size_t               index;
    typedef glutils::vector3f    point;
    typedef std::array<index, 3> index_list;
    typedef std::vector<point>   point_list;

    ///@}
    ///@name Internal State
    ///@{

    index_list m_indexes;       ///< The indexes in m_points held by this.
    const point_list& m_points; ///< A reference to the points referenced here.

    vector3f m_normal;          ///< The right-handed normal.

    ///@}

    public:

      ///@name Construction
      ///@{

      triangle_facet(const index _i1, const index _i2, const index _i3,
          const point_list& _pl) : m_indexes{_i1, _i2, _i3}, m_points(_pl) {
        compute_normal();
      }

      ///@}
      ///@name Accessors
      ///@{

      typedef index_list::const_iterator iterator;

      iterator begin() const noexcept {return m_indexes.begin();}
      iterator end() const noexcept {return m_indexes.end();}

      const index operator[](const size_t _i) const noexcept {
        return m_indexes[_i];
      }

      const point& get_point(const size_t _i) const noexcept {
        return m_points[_i];
      }

      const vector3f& get_normal() const noexcept {return m_normal;}

      ///@}
      ///@name Helpers
      ///@{

      void compute_normal() noexcept;

      ///@}

      friend class triangulated_model;

  };

  //////////////////////////////////////////////////////////////////////////////
  /// A geometric model, represented by vertices and triangular faces.
  //////////////////////////////////////////////////////////////////////////////
  class triangulated_model {

    ///@name Local Types
    ///@{

    typedef glutils::vector3f       point;
    typedef std::vector<point>      point_list;
    typedef glutils::triangle_facet facet;
    typedef std::vector<facet>      facet_list;

    ///@}
    ///@name Internal State
    ///@{

    point_list m_points;  ///< The vertices of this model.
    facet_list m_facets;  ///< The facets of this model.

    ///@}

    public:

      ///@name Creation Interface
      ///@{

      /// Add a point to the model.
      const size_t add_point(const point& _p, const bool _duplicates = false);

      /// Add a facet to the model by referencing the indexes of existing points.
      const size_t add_facet(const size_t _i1, const size_t _i2,
          const size_t _i3);

      ///@}
      ///@name Accessors
      ///@{

      const size_t num_points() const noexcept;
      const size_t num_facets() const noexcept;

      const point& get_point(const size_t _i) const noexcept;
      const facet& get_facet(const size_t _i) const noexcept;

      typedef point_list::const_iterator point_iterator;

      point_iterator points_begin() const noexcept {return m_points.begin();}
      point_iterator points_end() const noexcept {return m_points.end();}

      typedef facet_list::const_iterator facet_iterator;

      facet_iterator facets_begin() const noexcept {return m_facets.begin();}
      facet_iterator facets_end() const noexcept {return m_facets.end();}

      const vector3f find_center() const noexcept;

      ///@}
      ///@name Modifiers
      ///@{

      void translate(const vector3f& _v) noexcept;

      void clean() noexcept;

      ///@}
  };

}

#endif
