/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/intersection/intersection.hpp"

namespace detray::svgtools::utils {

/// @returns detray point of an intersection.
template <typename detector_t>
inline auto intersection_point(
    const typename detector_t::geometry_context& context,
    const detector_t& detector,
    const detray::intersection2D<typename detector_t::surface_type,
                                 typename detector_t::transform3>&
        d_intersection) {
    const typename detector_t::vector3 dir{};
    const detray::surface surface{detector, d_intersection.sf_desc};
    return surface.local_to_global(context, d_intersection.local, dir);
}

}  // namespace detray::svgtools::utils