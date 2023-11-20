/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/detector_volume.hpp"
#include "detray/plugins/svgtools/conversion/volume.hpp"

// Actsvg include(s)
#include "actsvg/proto/detector.hpp"

namespace detray::svgtools::conversion {

/// @brief Generates the proto detector object
///
/// @param context The geometry context.
/// @param detector The detector object.
/// @param hide_portals whether to display portals.
/// @param hide_passives whether to display passive surfaces.
///
/// @returns An actsvg proto volume representing the volume.
template <typename point3_container_t, typename detector_t, typename view_t>
auto detector(const typename detector_t::geometry_context& context,
              const detector_t& detector, const view_t& view,
              bool hide_portals = false, bool hide_passives = false,
              bool hide_grids = false) {

    actsvg::proto::detector<point3_container_t> p_detector;

    for (const auto& vol_desc : detector.volumes()) {

        p_detector._volumes.push_back(
            svgtools::conversion::volume<point3_container_t>(
                context, detector, detector_volume{detector, vol_desc}, view,
                hide_portals, hide_passives, hide_grids));
    }

    return p_detector;
}

}  // namespace detray::svgtools::conversion
