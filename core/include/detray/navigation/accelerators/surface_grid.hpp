/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/utils/grid/detail/axis_helpers.hpp"
#include "detray/utils/grid/detail/grid_bins.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/grid/grid_collection.hpp"
#include "detray/utils/type_traits.hpp"

namespace detray::detail {

// TODO: Define concepts

template <class accelerator_t>
struct is_surface_grid<
    accelerator_t,
    std::enable_if_t<
        is_grid<accelerator_t>::value  &&
            std::is_same<
                typename accelerator_t::value_type,
                surface_descriptor<
                    typename accelerator_t::value_type::mask_link,
                    typename accelerator_t::value_type::material_link,
                    typename accelerator_t::value_type::transform_link,
                    typename accelerator_t::value_type::navigation_link>>::value ,
        void>> : public std::true_type {};

}  // namespace detray::detail
