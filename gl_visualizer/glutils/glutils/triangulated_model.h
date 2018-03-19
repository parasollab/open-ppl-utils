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
  class triangle_facet final
  {

    ///@name Local Types
    ///@{

    typedef size_t                     index;
    typedef std::array<index, 3>       index_list;
    typedef index_list::const_iterator iterator;

    typedef glutils::vector3f          point;
    typedef std::vector<point>         point_list;

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
          const point_list& _pl);

      ///@}
      ///@name Accessors
      ///@{

      iterator begin() const noexcept; ///< Begin iterator over the facet indexes.
      iterator end() const noexcept;   ///< End iterator over the facet indexes.

      /// Get the index (in the owning model) of a facet point.
      /// @param[in] _i The index (in this facet) of the point.
      /// @return The index of point _i in the owning model.
      index operator[](const size_t _i) const noexcept;

      /// Get a facet point.
      /// @param[in] _i The index (in this facet) of the point.
      /// @return The point in the owning model referenced by facet index _i.
      const point& get_point(const size_t _i) const noexcept;

      /// Get the facet normal.
      const vector3f& get_normal() const noexcept;

      ///@}
      ///@name Helpers
      ///@{

      /// Compute the facet normal.
      void compute_normal() noexcept;

      ///@}
      ///@name Equality
      ///@{
      /// Determine whether two facets represent the same set of points in the
      /// same order.

      bool operator==(const triangle_facet& _t) const noexcept;
      bool operator!=(const triangle_facet& _t) const noexcept;

      ///@}
      ///@name Ordering
      ///@{

      /// Defines a weak ordering to allow sorting.
      bool operator<(const triangle_facet& _t) const noexcept;

      ///@}

      friend class triangulated_model;

  };

  //////////////////////////////////////////////////////////////////////////////
  /// A geometric model, represented by vertices and triangular faces.
  //////////////////////////////////////////////////////////////////////////////
  class triangulated_model final
  {

    ///@name Local Types
    ///@{

    typedef glutils::vector3f          point;
    typedef std::vector<point>         point_list;
    typedef point_list::const_iterator point_iterator;

    typedef glutils::triangle_facet    facet;
    typedef std::vector<facet>         facet_list;
    typedef facet_list::const_iterator facet_iterator;

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
      /// @param[in] _p The point to add.
      /// @param[in] _duplicates Allow duplicate points?
      /// @return The index of the added point.
      size_t add_point(const point& _p, const bool _duplicates = false);

      /// Add a facet to the model by referencing the indexes of existing points.
      /// @param[in] _i1 The index of the first point.
      /// @param[in] _i1 The index of the second point.
      /// @param[in] _i1 The index of the third point.
      /// @return The index of the added facet.
      size_t add_facet(const size_t _i1, const size_t _i2, const size_t _i3);

      ///@}
      ///@name Accessors
      ///@{

      const point& get_point(const size_t _i) const noexcept;
      const facet& get_facet(const size_t _i) const noexcept;

      point_iterator points_begin() const noexcept;
      point_iterator points_end() const noexcept;

      facet_iterator facets_begin() const noexcept;
      facet_iterator facets_end() const noexcept;

      ///@}
      ///@name Queries
      ///@{

      /// Get the number of vertices in the model.
      size_t num_points() const noexcept;

      /// Get the number of facets in the model.
      size_t num_facets() const noexcept;

      /// Find the model's centroid by averaging all of its vertices.
      /// @return The centroid of the model's vertices.
      const vector3f find_center() const noexcept;

      ///@}
      ///@name Modifiers
      ///@{

      /// Translate all of the vertices in the model.
      /// @param[in] _v The translation vector to apply.
      void translate(const vector3f& _v) noexcept;

      /// Remove duplicate vertices from the model.
      void clean() noexcept;

      ///@}
      ///@name Equality
      ///@{
      /// Determine whether two models contain the same facets and vertices.

      bool operator==(const triangulated_model& _t) const noexcept;
      bool operator!=(const triangulated_model& _t) const noexcept;

      ///@}
      ///@name Common Shapes
      ///@{

      /// Create a box model.
      /// @param[in] _lenX The box length in the x direction.
      /// @param[in] _lenY The box length in the y direction.
      /// @param[in] _lenZ The box length in the z direction.
      /// @return A triangulated box model of size _lenX x _lenY x _lenZ.
      static triangulated_model make_box(GLfloat _lenX = 1, GLfloat _lenY = 1,
          GLfloat _lenZ = 1);

      /// Create a sphere centered at the origin with axis along the z direction.
      /// @param[in] _radius The sphere radius.
      /// @param[in] _segments The number of segments to use. The model will
      ///                      have _segments rings of _segments squares, or
      ///                      2 * _segments^2 triangles.
      /// @return A triangulated sphere model of radius _radius.
      static triangulated_model make_sphere(const GLfloat _radius = 1,
          const size_t _segments = 16);

      /// Draw a cone with the base centered at the origin and tip pointed away
      /// from the camera.
      /// @param[in] _radius The radius of the base.
      /// @param[in] _height The height of the cone.
      /// @param[in] _segments The number of segments to use for the sides.
      static triangulated_model make_cone(const GLfloat _radius = 1,
          const GLfloat _height = 1, const size_t _segments = 16);

      /// Draw a cylinder centered at the origin and oriented along the z-axis.
      /// @param[in] _radius The cylinder radius.
      /// @param[in] _length The length perpendicular to the radius.
      /// @param[in] _segments The number of segments to use for the side wall.
      static triangulated_model make_cylinder(const GLfloat _radius = 1,
          const GLfloat _length = 1, const size_t _segments = 16);

      ///@}

  };

}

#endif
