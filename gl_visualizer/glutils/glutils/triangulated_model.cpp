#include "glutils/triangulated_model.h"

#include <algorithm>
#include <list>
#include <map>

namespace glutils {

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Triangle Facet ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /*----------------------------- Construction -------------------------------*/

  triangle_facet::
  triangle_facet(const index _i1, const index _i2, const index _i3,
      const point_list& _pl)
    : m_indexes{_i1, _i2, _i3}, m_points(_pl)
  {
    // Rotate the vertex list so that the lowest index is always first. This
    // helps check for equality efficiently.
    auto iter = std::min_element(m_indexes.begin(), m_indexes.end());
    std::rotate(m_indexes.begin(), iter, m_indexes.end());
    compute_normal();
  }

  /*----------------------------- Accessors ----------------------------------*/

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

  /*------------------------------ Helpers -----------------------------------*/

  void
  triangle_facet::
  compute_normal() noexcept
  {
    m_normal = (get_point(1) - get_point(0)) % (get_point(2) - get_point(0));
    m_normal.normalize();
  }

  /*------------------------------ Equality ----------------------------------*/

  bool
  triangle_facet::
  operator==(const triangle_facet& _t) const noexcept
  {
    // If both facets have the same point list, we can just check the indexes.
    if(&m_points == &_t.m_points) {
      for(size_t i = 0; i < 3; ++i)
        if(m_indexes[i] != _t.m_indexes[i])
          return false;
    }
    // Otherwise, we need to check the actual points.
    else {
      for(size_t i = 0; i < 3; ++i)
        if(get_point(i) != _t.get_point(i))
          return false;
    }
    return true;
  }


  bool
  triangle_facet::
  operator!=(const triangle_facet& _t) const noexcept
  {
    return !(*this == _t);
  }

  /*------------------------------ Ordering ----------------------------------*/

  bool
  triangle_facet::
  operator<(const triangle_facet& _t) const noexcept
  {
    for(size_t i = 0; i < 3; ++i) {
      if(m_indexes[0] < _t.m_indexes[0])
        return true;
      else if(_t.m_indexes[0] < m_indexes[0])
        return false;
    }
    return false;
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
    /// @TODO Fix this silly O(n^2) algorithm and replace with O(n) solution by
    ///       re-writing a new point list and swapping.
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

  /*------------------------------ Equality ----------------------------------*/

  bool
  triangulated_model::
  operator==(const triangulated_model& _t) const noexcept
  {
    // First check that the number of points are the same.
    if(num_points() != _t.num_points())
      return false;

    // If the number of points are the same, test facets.
    // ARG facets can be in a different order @#$%
    // ARG referenced points can be in a different order @#$%
    return m_facets == _t.m_facets;
  }


  bool
  triangulated_model::
  operator!=(const triangulated_model& _t) const noexcept
  {
    return !(*this == _t);
  }

  /*---------------------------- Common Shapes -------------------------------*/

  triangulated_model
  triangulated_model::
  make_box(GLfloat _lenX, GLfloat _lenY, GLfloat _lenZ)
  {
    // Get half-lengths.
    _lenX /= 2;
    _lenY /= 2;
    _lenZ /= 2;

    // Make vertices.
    vector3f pts[8] = {{-_lenX,  _lenY,  _lenZ},
                       {-_lenX, -_lenY,  _lenZ},
                       { _lenX, -_lenY,  _lenZ},
                       { _lenX,  _lenY,  _lenZ},
                       {-_lenX,  _lenY, -_lenZ},
                       {-_lenX, -_lenY, -_lenZ},
                       { _lenX, -_lenY, -_lenZ},
                       { _lenX,  _lenY, -_lenZ}};

    // Add points to the model and get their indexes.
    triangulated_model t;
    size_t id[8];
    for(size_t i = 0; i < 8; ++i)
      id[i] = t.add_point(pts[i]);

    // Make facets.
    size_t facets[12][3] = {{id[0], id[1], id[2]},
                            {id[0], id[2], id[3]},
                            {id[3], id[2], id[6]},
                            {id[3], id[6], id[7]},
                            {id[7], id[6], id[5]},
                            {id[7], id[5], id[4]},
                            {id[4], id[5], id[1]},
                            {id[4], id[1], id[0]},
                            {id[0], id[3], id[7]},
                            {id[0], id[7], id[4]},
                            {id[1], id[5], id[6]},
                            {id[1], id[6], id[2]}};

    // Add facets.
    for(size_t i = 0; i < 12; ++i)
      t.add_facet(facets[i][0], facets[i][1], facets[i][2]);

    return t;
  }


  triangulated_model
  triangulated_model::
  make_sphere(const GLfloat _radius, const size_t _segments)
  {
    triangulated_model t;

    const GLfloat zIncr = glutils::PI / _segments; // Angle increment for z.
    const GLfloat oIncr = 2 * zIncr;               // Angle increment for x,y.

    GLfloat x, y, z, r;

    // Draw +zHat cap.
    {
      // Create the +zHat pole.
      x = 0;
      y = 0;
      z = _radius;
      const size_t capIndex = t.add_point({x, y, z});

      // Create the ring of points beneath the pole.
      std::vector<size_t> indexes;;

      z = _radius * std::cos(zIncr);
      r = _radius * std::sin(zIncr);
      for(size_t i = 0; i < _segments; ++i) {
        x = r * std::cos(oIncr * i);
        y = r * std::sin(oIncr * i);
        indexes.push_back(t.add_point({x, y, z}));
      }

      // Create the facets connecting the pole to the point ring.
      for(size_t i = 0; i < _segments; ++i)
        t.add_facet(capIndex, indexes[i], indexes[(i + 1) % _segments]);
    }

    // Draw main surface.
    {
      GLfloat z2, r2;

      // Create a ring of segments following the previous.
      std::vector<size_t> topIndexes, bottomIndexes;

      for(size_t j = 1; j < _segments - 1; ++j) {
        // The top and bottom point rings in this segment ring each require a
        // different z and planar radius.
        z  = _radius * std::cos(zIncr * j);
        r  = _radius * std::sin(zIncr * j);
        z2 = _radius * std::cos(zIncr * (j + 1));
        r2 = _radius * std::sin(zIncr * (j + 1));

        // Generate the points for this segment ring.
        topIndexes.clear();
        bottomIndexes.clear();
        for(size_t i = 0; i < _segments; ++i) {
          x = std::cos(oIncr * i);
          y = std::sin(oIncr * i);
          topIndexes.push_back(   t.add_point({x * r , y * r ,  z}));
          bottomIndexes.push_back(t.add_point({x * r2, y * r2, z2}));
        }

        // Create facets to complete this segment ring.
        for(size_t i = 0; i < _segments; ++i) {
          const size_t i2 = (i + 1) % _segments;
          t.add_facet(topIndexes[i] , bottomIndexes[i], topIndexes[i2]);
          t.add_facet(topIndexes[i2], bottomIndexes[i], bottomIndexes[i2]);
        }
      }
    }

    // Draw -zHat cap.
    {
      // Create the +zHat pole.
      x = 0;
      y = 0;
      z = -_radius;
      const size_t capIndex = t.add_point({x, y, z});

      // Create the ring of points above the pole.
      std::vector<size_t> indexes;

      z = _radius * std::cos(zIncr * (_segments - 1));
      r = _radius * std::sin(zIncr * (_segments - 1));
      for(size_t i = 0; i < _segments; ++i) {
        x = r * std::cos(oIncr * i);
        y = r * std::sin(oIncr * i);
        indexes.push_back(t.add_point({x, y, z}));
      }

      // Create the facets connecting the pole to the point ring.
      for(size_t i = 0; i < _segments; ++i)
        t.add_facet(capIndex, indexes[(i + 1) % _segments], indexes[i]);
    }

    return t;
  }


  triangulated_model
  triangulated_model::
  make_cone(const GLfloat _radius, const GLfloat _height, const size_t _segments)
  {
    triangulated_model t;

    const GLfloat incr = 2. * glutils::PI / _segments;
    GLfloat x, y, z;

    // Create the point ring.
    z = 0;
    std::vector<size_t> ringIndexes;
    for(size_t i = 0; i < _segments; ++i) {
      x = _radius * std::cos(incr * i);
      y = _radius * std::sin(incr * i);
      ringIndexes.push_back(t.add_point({x, y, z}));
    }

    // Create the circular base on the x-y plane.
    {
      // Create the cap.
      x = 0;
      y = 0;
      z = 0;
      const size_t capIndex = t.add_point({x, y, z});

      // Create facets connecting the cap to the point ring.
      for(size_t i = 0; i < _segments; ++i)
        t.add_facet(capIndex, ringIndexes[i], ringIndexes[(i + 1) % _segments]);
    }

    // Create the conic cap on top of the circle.
    {
      // Create the cap.
      x = 0;
      y = 0;
      z = -_height;
      const size_t capIndex = t.add_point({x, y, z});

      // Create facets connecting the cap to the point ring.
      for(size_t i = 0; i < _segments; ++i)
        t.add_facet(capIndex, ringIndexes[(i + 1) % _segments], ringIndexes[i]);
    }

    // Translate the model so that it is centered on its bounding box.
    t.translate({0, 0, _height/2});
    return t;
  }


  triangulated_model
  triangulated_model::
  make_cylinder(const GLfloat _radius, const GLfloat _length,
      const size_t _segments)
  {
    triangulated_model t;

    const GLfloat incr = 2. * glutils::PI / _segments,
                  halfLength = _length / 2.;
    GLfloat x, y, z;

    // Draw top cap.
    x = 0;
    y = 0;
    z = halfLength;
    const size_t topCapIndex = t.add_point({x, y, z});

    // Create the top point ring.
    std::vector<size_t> topIndexes;
    for(size_t i = 0; i < _segments; ++i) {
      x = _radius * std::cos(incr * i);
      y = _radius * std::sin(incr * i);
      topIndexes.push_back(t.add_point({x, y, z}));
    }

    // Create facets connecting the cap to the point ring.
    for(size_t i = 0; i < _segments; ++i)
      t.add_facet(topCapIndex, topIndexes[i], topIndexes[(i + 1) % _segments]);

    // Draw bottom cap.
    x = 0;
    y = 0;
    z = -halfLength;
    const size_t bottomCapIndex = t.add_point({x, y, z});

    // Create the bottom point ring.
    std::vector<size_t> bottomIndexes;
    for(size_t i = 0; i < _segments; ++i) {
      x = _radius * std::cos(incr * i);
      y = _radius * std::sin(incr * i);
      bottomIndexes.push_back(t.add_point({x, y, z}));
    }

    // Create facets connecting the cap to the point ring.
    for(size_t i = 0; i < _segments; ++i)
      t.add_facet(bottomCapIndex, bottomIndexes[(i + 1) % _segments],
                  bottomIndexes[i]);

    // Create facets connecting the top and bottom rings.
    for(size_t i = 0; i < _segments; ++i) {
      const size_t i2 = (i + 1) % _segments;
      t.add_facet(topIndexes[i] , bottomIndexes[i], topIndexes[i2]);
      t.add_facet(topIndexes[i2], bottomIndexes[i], bottomIndexes[i2]);
    }

    return t;
  }

  /*--------------------------------------------------------------------------*/
}
