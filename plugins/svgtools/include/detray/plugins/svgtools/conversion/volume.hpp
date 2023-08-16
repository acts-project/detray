/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/detector_volume.hpp"
#include "detray/plugins/svgtools/conversion/portal.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/utils/volume_utils.hpp"

// Actsvg include(s)
#include "actsvg/proto/portal.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools::conversion {

/// @brief Calculates the proto volume of a collection of detray surfaces.
///
/// @param context The context.
/// @param d_surfaces The detray surfaces.
///
/// @returns An actsvg proto volume representing the volume.
template <typename point3_container_t, typename detector_t>
auto volume(const typename detector_t::geometry_context& context,
            const detector_t& detector,
            const std::vector<detray::surface<detector_t>>& d_surfaces) {
    using p_volume_t = actsvg::proto::volume<point3_container_t>;
    using p_portal_t = actsvg::proto::portal<point3_container_t>;
    using p_surface_t = actsvg::proto::surface<point3_container_t>;
    p_volume_t p_volume;
    std::vector<p_surface_t> surfaces;
    std::vector<p_portal_t> portals;
    for (const auto& item : d_surfaces) {
        if (item.is_portal()) {
            auto portal = svgtools::conversion::portal<point3_container_t>(
                context, detector, item);
            portals.push_back(portal);
        } else {
            auto surface = svgtools::conversion::surface<point3_container_t>(
                context, item);
            surfaces.push_back(surface);
        }
    }
    p_volume._v_surfaces = surfaces;
    p_volume._portals = portals;
    return p_volume;
}

/// @brief Calculates the proto volume of a detray volume.
///
/// @param context The context.
/// @param d_volume The detray volume.
///
/// @returns An actsvg proto volume representing the volume.
template <typename point3_container_t, typename detector_t>
auto volume(const typename detector_t::geometry_context& context,
            const detector_t& detector,
            const detray::detector_volume<detector_t>& d_volume) {
    std::vector<detray::surface<detector_t>> surfaces{};
    const auto descriptors =
        svgtools::utils::surface_lookup(detector, d_volume);
    for (const auto& desc : descriptors) {
        surfaces.push_back(detray::surface{detector, desc});
    }
    return svgtools::conversion::volume<point3_container_t>(context, detector,
                                                            surfaces);
}

}  // namespace detray::svgtools::conversion