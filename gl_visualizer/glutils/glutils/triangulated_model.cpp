#include "glutils/triangulated_model.h"

#include <algorithm>
#include <list>
#include <map>

namespace glutils {

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Triangle Facet ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  triangle_facet::
  triangle_facet(const index _i1, const index _i2, const index _i3,
      const point_list& _pl)
    : m_indexes{_i1, _i2, _i3}, m_points(_pl)
  {
    compute_normal();
  }


  triangle_facet::iterator
  triangle_facet::
  begin() const noexcept
  {
    return m_indexes.begin();
  }


  triangle_facet::iterator
  triangle_facet::
  end() const noexcept
  {
    return m_indexes.end();
  }


  triangle_facet::index
  triangle_facet::
  operator[](const size_t _i) const noexcept
  {
    return m_indexes[_i];
  }


  const triangle_facet::point&
  triangle_facet::
  get_point(const size_t _i) const noexcept
  {
    return m_points[_i];
  }


  const vector3f&
  triangle_facet::
  get_normal() const noexcept
  {
    return m_normal;
  }


  void
  triangle_facet::
  compute_normal() noexcept
  {
    m_normal = (get_point(1) - get_point(0)) % (get_point(2) - get_point(0));
    m_normal.normalize();
  }

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~ Triangulated Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /*--------------------------- Creation Interface ---------------------------*/

  size_t
  triangulated_model::
  add_point(const point& _p, const bool _duplicates)
  {
    if(!_duplicates) {
      auto iter = std::find(m_points.begin(), m_points.end(), _p);
      if(iter != m_points.end())
        return std::distance(m_points.begin(), iter);
    }
    m_points.push_back(_p);
    return m_points.size() - 1;
  }


  size_t
  triangulated_model::
  add_facet(const size_t _i1, const size_t _i2, const size_t _i3)
  {
    m_facets.emplace_back(_i1, _i2, _i3, m_points);
    return m_facets.size() - 1;
  }

  /*------------------------------- Accessors --------------------------------*/

  const triangulated_model::point&
  triangulated_model::
  get_point(const size_t _i) const noexcept
  {
    return m_points[_i];
  }


  const triangulated_model::facet&
  triangulated_model::
  get_facet(const size_t _i) const noexcept
  {
    return m_facets[_i];
  }


  triangulated_model::point_iterator
  triangulated_model::
  points_begin() const noexcept
  {
    return m_points.begin();
  }


  triangulated_model::point_iterator
  triangulated_model::
  points_end() const noexcept
  {
    return m_points.end();
  }


  triangulated_model::facet_iterator
  triangulated_model::
  facets_begin() const noexcept
  {
    return m_facets.begin();
  }


  triangulated_model::facet_iterator
  triangulated_model::
  facets_end() const noexcept
  {
    return m_facets.end();
  }

  /*------------------------------- Queries ----------------------------------*/

  size_t
  triangulated_model::
  num_points() const noexcept
  {
    return m_points.size();
  }


  size_t
  triangulated_model::
  num_facets() const noexcept
  {
    return m_facets.size();
  }


  const vector3f
  triangulated_model::
  find_center() const noexcept
  {
    vector3f center;
    for(const auto& p : m_points)
      center += p;
    return center /= m_points.size();
  }

  /*------------------------------ Modifiers ---------------------------------*/

  void
  triangulated_model::
  translate(const vector3f& _v) noexcept
  {
    for(auto& p : m_points)
      p += _v;
  }


  void
  triangulated_model::
  clean() noexcept
  {
    // Remove vertices that don't refer to anything.

    // Create a map from each vertex index to the facets that contain it.
    using std::map;
    using std::list;
    map<size_t, list<facet*>> index_map;

    // Add each index to the map with an empty facet list.
    for(size_t i = 0; i < num_points(); ++i)
      index_map[i];

    // Iterate through facets and add them to the map.
    for(auto iter = m_facets.begin(); iter != m_facets.end(); ++iter)
      for(const auto i : *iter)
        index_map[i].push_back(&*iter);

    // Remove unused vertices, starting from the back to avoid having to
    // recompute the indexes.
    for(auto iter = index_map.rbegin(); iter != index_map.rend(); ++iter) {
      const bool unused = iter->second.empty();
      if(unused)
        m_points.erase(m_points.begin() + iter->first);
    }

    // When we removed the unused vertices, we shifted the indexes of the
    // remaining points. We will need to adjust the indexes in all of the facets
    // appropriately. First, figure out how much each index needs to change.
    map<size_t, size_t> change_map;
    size_t offset = 0;
    for(const auto& pair : index_map) {
      const bool unused = pair.second.empty();
      offset += unused;
      if(!unused)
        change_map[pair.first] = offset;
    }

    // Now iterate through the facets and adjust the indexes using the change
    // map.
    for(auto iter = m_facets.begin(); iter != m_facets.end(); ++iter)
      for(auto& index : iter->m_indexes)
        index -= change_map.at(index);
  }

  /*--------------------------------------------------------------------------*/
}
