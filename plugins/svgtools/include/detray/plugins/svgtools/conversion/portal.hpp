/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/svgtools/conversion/link.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/utils/link_utils.hpp"

// Actsvg includes(s)
#include "actsvg/proto/portal.hpp"

namespace detray::svgtools::conversion {

/// @returns An actsvg proto portal representing the portal.
/// @note detray portal is_portal() should be true.
template <typename point3_container_t, typename detector_t>
auto portal(const typename detector_t::geometry_context& context,
            const detector_t& detector,
            const detray::surface<detector_t>& d_portal) {
    assert(d_portal.is_portal());
    using p_portal_t = actsvg::proto::portal<point3_container_t>;
    p_portal_t p_portal;
    if (svgtools::utils::is_not_world_portal(d_portal)) {
        p_portal._volume_links = {
            svgtools::conversion::link<point3_container_t>(context, detector,
                                                           d_portal)};
    }
    p_portal._surface =
        svgtools::conversion::surface<point3_container_t>(context, d_portal);
    return p_portal;
}

}  // namespace detray::svgtools::conversion