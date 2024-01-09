/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/detail/grid_bins.hpp"
#include "detray/surface_finders/grid/grid.hpp"
#include "detray/surface_finders/grid/grid_collection.hpp"
#include "detray/utils/type_traits.hpp"

namespace detray::detail {

template <typename multi_axis_t, typename bin_t,
          template <std::size_t> class serializer_t>
struct is_grid<grid<multi_axis_t, bin_t, serializer_t>>
    : public std::true_type {};

}  // namespace detray::detail
