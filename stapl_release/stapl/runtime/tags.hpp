/*
// Copyright (c) 2000-2009, Texas Engineering Experiment Station (TEES), a
// component of the Texas A&M University System.

// All rights reserved.

// The information and source code contained herein is the exclusive
// property of TEES and may not be disclosed, examined or reproduced
// in whole or in part without explicit written authorization from TEES.
*/


#ifndef STAPL_RUNTIME_TAGS_HPP
#define STAPL_RUNTIME_TAGS_HPP

#include "config/types.hpp"
#include <limits>

namespace stapl {

//////////////////////////////////////////////////////////////////////
/// @brief Tag type to represent all locations of a gang.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
struct all_locations_t
{ };

//////////////////////////////////////////////////////////////////////
/// @brief Tag to represent all locations of a gang.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
constexpr all_locations_t all_locations = { };


//////////////////////////////////////////////////////////////////////
/// @brief Tag type to represent the current location.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
struct this_location_t
{ };

//////////////////////////////////////////////////////////////////////
/// @brief Tag to represent the current location in a gang.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
constexpr this_location_t this_location = { };


//////////////////////////////////////////////////////////////////////
/// @brief Tag type to represent the root of a collective operation.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
struct root_location_t
{ };

//////////////////////////////////////////////////////////////////////
/// @brief Tag to represent the root of a collective operation.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
constexpr root_location_t root_location = { };

//////////////////////////////////////////////////////////////////////
/// @brief Tag type to represent the neighbor locations of the current location.
///
/// Neighbor locations are the ones that are on the same hierarchy level as the
/// current.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
struct neighbor_locations_t
{ };

//////////////////////////////////////////////////////////////////////
/// @brief Tag to represent the neighbor locations of the current location.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
constexpr neighbor_locations_t neighbor_locations = { };


//////////////////////////////////////////////////////////////////////
/// @brief Tag type for an invalid argument.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
struct none_t
{ };

//////////////////////////////////////////////////////////////////////
/// @brief Tag for invalid argument.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
constexpr none_t none = { };


//////////////////////////////////////////////////////////////////////
/// @brief Tag type to execute all scheduled entries.
///
/// @ingroup executors
//////////////////////////////////////////////////////////////////////
struct execute_all_t
{ };

//////////////////////////////////////////////////////////////////////
/// @brief Tag to execute all scheduled entries.
///
/// @ingroup executors
//////////////////////////////////////////////////////////////////////
constexpr execute_all_t execute_all = { };


//////////////////////////////////////////////////////////////////////
/// @brief Tag type to describe which level of the hierarchy is requested.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
class level
{
public:
  static const level_type invalid = std::numeric_limits<level_type>::max();
  static const level_type max     = std::numeric_limits<level_type>::max() - 1;
  static const level_type current = std::numeric_limits<level_type>::max() - 2;
  static const level_type lower   = std::numeric_limits<level_type>::max() - 3;

private:
  level_type m_level;
  bool       m_restricted;

public:
  constexpr level(void)
  : m_level(invalid),
    m_restricted(false)
  { }

  constexpr level(this_location_t)
  : m_level(max),
    m_restricted(true)
  { }

  constexpr level(all_locations_t)
  : m_level(current),
    m_restricted(false)
  { }

  constexpr level(const level_type lvl, const bool restricted = false)
  : m_level(lvl),
    m_restricted(restricted)
  { }

  constexpr level_type get_level(void) const noexcept
  { return m_level; }

  constexpr bool is_restricted(void) const noexcept
  { return m_restricted; }
};


constexpr bool operator==(level const& x, level const& y)
{
  return (x.get_level()==y.get_level() && x.is_restricted()==y.is_restricted());
}

constexpr bool operator!=(level const& x, level const& y)
{
  return !(x==y);
}


//////////////////////////////////////////////////////////////////////
/// @brief Tag for lowest level.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
constexpr level lowest_level{level::max, true};

//////////////////////////////////////////////////////////////////////
/// @brief Tag for current level.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
constexpr level current_level{level::current, false};

//////////////////////////////////////////////////////////////////////
/// @brief Tag for one level down.
///
/// @ingroup ARMITags
//////////////////////////////////////////////////////////////////////
constexpr level lower_level{level::lower, true};


namespace runtime {

//////////////////////////////////////////////////////////////////////
/// @brief Tag type to defer an operation.
///
/// @ingroup runtimeMetadata
//////////////////////////////////////////////////////////////////////
struct deferred_t
{ };

//////////////////////////////////////////////////////////////////////
/// @brief Tag to create @ref gang_md objects with deferred id.
///
/// @ingroup runtimeMetadata
//////////////////////////////////////////////////////////////////////
constexpr deferred_t deferred = { };


//////////////////////////////////////////////////////////////////////
/// @brief Tag type to skip any operations related to a @ref context.
///
/// @ingroup runtimeMetadata
//////////////////////////////////////////////////////////////////////
struct no_context_t
{ };

//////////////////////////////////////////////////////////////////////
/// @brief Tag to skip any operations related to a @ref context.
///
/// @ingroup runtimeMetadata
//////////////////////////////////////////////////////////////////////
constexpr no_context_t no_context = { };


////////////////////////////////////////////////////////////////////
/// @brief Tag type to mark that a flush is not implicit.
///
/// @ingroup requestBuildingBlock
////////////////////////////////////////////////////////////////////
struct no_implicit_flush_t
{ };

////////////////////////////////////////////////////////////////////
/// @brief Tag to mark that a flush is not implicit.
///
/// @ingroup requestBuildingBlock
////////////////////////////////////////////////////////////////////
constexpr no_implicit_flush_t no_implicit_flush = { };

} // namespace runtime

} // namespace stapl

#endif

