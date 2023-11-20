/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/svgtools/conversion/grid.hpp"
#include "detray/plugins/svgtools/conversion/portal.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/utils/volume_utils.hpp"

// Actsvg include(s)
#include "actsvg/proto/volume.hpp"

namespace detray::svgtools::conversion {

/// @brief Generates the proto volume of a detray volume.
///
/// @param context The geometry context.
/// @param detector The detractor the volume belongs to.
/// @param d_volume The detray volume.
/// @param hide_portals whether to display the volumes portals.
/// @param hide_passives whether to display the contained passive surfaces.
///
/// @returns An actsvg proto volume representing the volume.
template <typename point3_container_t, typename detector_t, typename view_t>
auto volume(const typename detector_t::geometry_context& context,
            const detector_t& detector,
            const detray::detector_volume<detector_t>& d_volume,
            const view_t& view, bool hide_portals = false,
            bool hide_passives = false, bool hide_grids = false) {

    actsvg::proto::volume<point3_container_t> p_volume;
    p_volume._index = d_volume.index();

    for (const auto& desc :
         svgtools::utils::surface_lookup(detector, d_volume)) {

        const auto sf = detray::surface{detector, desc};

        if (sf.is_portal()) {
            if (!hide_portals) {
                auto portal = svgtools::conversion::portal<point3_container_t>(
                    context, detector, sf, false);
                p_volume._portals.push_back(portal);
            }
        } else if (!(sf.is_passive() && hide_passives)) {
            auto surface =
                svgtools::conversion::surface<point3_container_t>(context, sf);
            p_volume._v_surfaces.push_back(surface);
        }
    }

    // Convert grid, if present
    if (!hide_grids) {
        if (auto p_grid_ptr = svgtools::conversion::grid<actsvg::scalar>(
                detector, p_volume._index, view)) {
            p_volume._surface_grid = *p_grid_ptr;
        }
    }

    return p_volume;
}

}  // namespace detray::svgtools::conversion
