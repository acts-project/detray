/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/grid.hpp"
#include "detray/surface_finders/grid/grid_collection.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/surface_finders/grid/serializer.hpp"
#include "detray/utils/type_traits.hpp"

namespace detray::detail {

template <typename multi_axis_t, typename value_t,
          template <std::size_t> class serializer_t, typename populator_impl_t>
struct is_grid<grid<multi_axis_t, value_t, serializer_t, populator_impl_t>>
    : public std::true_type {};

}  // namespace detray::detail
