/*
// Copyright (c) 2000-2009, Texas Engineering Experiment Station (TEES), a
// component of the Texas A&M University System.

// All rights reserved.

// The information and source code contained herein is the exclusive
// property of TEES and may not be disclosed, examined or reproduced
// in whole or in part without explicit written authorization from TEES.
*/

#ifndef STAPL_CONTAINERS_UNORDERED_SET_TRAITS_HPP
#define STAPL_CONTAINERS_UNORDERED_SET_TRAITS_HPP

#include <stapl/containers/distribution/directory/unordered_map_manager.hpp>
#include <stapl/containers/distribution/directory/unordered_registry.hpp>
#include <stapl/containers/distribution/directory/vector_directory.hpp>
#include <stapl/containers/distribution/container_manager/unordered_associative_container_manager.hpp>

#include "unordered_set_fwd.hpp"
#include <stapl/containers/unordered_set/distribution.hpp>
#include <stapl/containers/unordered_set/base_container.hpp>
#include <stapl/containers/unordered_set/base_container_traits.hpp>

#include <stapl/containers/unordered_map/unordered_map_partition.hpp>
#include <stapl/containers/mapping/mapper.hpp>
#include <stapl/domains/continuous.hpp>
#include <stapl/domains/iterator_domain.hpp>

namespace stapl {

//////////////////////////////////////////////////////////////////////
/// @brief Default traits for the unordered_set container. Specifies
///   customizable type parameters that could be changed
///   on a per-container basis.
/// @ingroup punorderedsetTraits
/// @tparam Key Type of the key. Each element has a unique key value
/// @tparam Hash Unary Functor returning a unique value of type size_t. The
///   only argument is a type Key.
/// @tparam Pred Binary predicate used to compare key values of elements stored
///   in the container for equality.
/// @tparam PS Partitioner type. Used to specify how to partition the elements.
/// @tparam M Mapper type. Maps the different partitions to the available
///   locations.
/// @see unordered_set
////////////////////////////////////////////////////////////////////////
template<typename Key, typename Hash, typename Pred, typename PS, typename M>
struct unordered_set_traits
{
  typedef Key                                                key_type;
  typedef Key                                                mapped_type;
  typedef Key                                                value_type;
  typedef Key                                                gid_type;
  typedef Hash                                               hasher;
  typedef Pred                                               key_equal;
  typedef PS                                                 partition_type;
  typedef M                                                  mapper_type;

  typedef unordered_set_base_container_traits<
    Key, Hash, Pred>                                 base_container_traits_type;
  typedef unordered_set_base_container<
    Key, base_container_traits_type>                 base_container_type;
  typedef unordered_map_manager<
    partition_type, mapper_type>                     manager_type;
  typedef vector_directory<
    partition_type, mapper_type, manager_type,
    unordered_registry<manager_type> >               directory_type;
  typedef unordered_associative_container_manager<
    base_container_type, manager_type>               container_manager_type;

  using container_type
    = unordered_set<Key, Hash, Pred, PS, M, unordered_set_traits>;
  using distribution_type = unordered_set_distribution<container_type>;
  using domain_type =
    iterator_domain<distribution_type, domain_impl::f_deref<gid_type>>;

  //////////////////////////////////////////////////////////////////////
  /// @brief Metafunction to compute the distribution type based on a
  /// container type.
  /// @tparam C Type of the container.
  //////////////////////////////////////////////////////////////////////
  template <typename C>
  struct construct_distribution
  {
    typedef unordered_set_distribution<C>  type;
  };
}; // struct unordered_set_traits

namespace unordered_set_impl {

//////////////////////////////////////////////////////////////////////
/// @brief Metafunction to compute traits class for the @ref unordered
///   set container. For primary template, some optional parameters are
///   defined, but no customized traits.  Use @ref unordered_set_traits.
/// @ingroup punorderedsetTraits
//////////////////////////////////////////////////////////////////////
template<typename Key, typename ...OptionalParams>
struct compute_unordered_set_traits
{
private:
  typedef tuple<
    stapl::hash<Key>,            // Hasher
    stapl::equal_to<Key>,        // Predicate
    unordered_map_partition<
      stapl::hash<Key>,
      continuous_domain<Key> >,  // Partition
    mapper<size_t>               // Mapper
  >                                                       default_value_types;

  typedef typename compute_type_parameters<
    default_value_types, OptionalParams...
  >::type                                                 param_types;

  typedef typename tuple_element<0, param_types>::type    hash_t;
  typedef typename tuple_element<1, param_types>::type    pred_t;
  typedef typename tuple_element<2, param_types>::type    partition_t;
  typedef typename tuple_element<3, param_types>::type    mapper_t;

public:
  typedef unordered_set_traits<Key, hash_t, pred_t, partition_t, mapper_t> type;
};


//////////////////////////////////////////////////////////////////////
/// @brief Specialization when traits is actually passed at class
//    template instantiation. Use passed type.
/// @ingroup punorderedsetTraits
//////////////////////////////////////////////////////////////////////
template<typename Key, typename Hash, typename Pred, typename PS, typename M,
  typename Traits>
struct compute_unordered_set_traits<Key, Hash, Pred, PS, M, Traits>
{
  typedef Traits type;
};


//////////////////////////////////////////////////////////////////////
/// @brief Specialization when none of the optional parameters are
///   provided at template instantiation.  Use @ref unordered_set_traits, with
///   default partition and mapper argument types.
/// @ingroup punorderedsetTraits
//////////////////////////////////////////////////////////////////////
template<typename Key>
struct compute_unordered_set_traits<Key>
{
  typedef unordered_set_traits<
    Key, stapl::hash<Key>, stapl::equal_to<Key>,
    unordered_map_partition<stapl::hash<Key>, continuous_domain<Key> >,
    mapper<size_t>
  > type;
};

} // namespace unordered_set_impl

} // namespace stapl

#endif
