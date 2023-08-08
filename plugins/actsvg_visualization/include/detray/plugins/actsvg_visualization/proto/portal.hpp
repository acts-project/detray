#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/proto/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/link.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/link_utils.hpp"

namespace detray::actsvg_visualization::proto {

/// @returns An actsvg proto portal representing the portal.
/// @note detray portal is_portal() must be true.
template <typename detector_t>
auto portal(const typename detector_t::geometry_context& context, const detector_t& detector, const detray::surface<detector_t>& d_portal)
{
    assert(d_portal.is_portal());
    proto_portal p_portal;
    if (utils::has_link(d_portal))
    {
        p_portal._volume_links = {proto::link(context, detector, d_portal)};
    }
    p_portal._surface = surface(context, d_portal);
    return p_portal;
}

}  // namespace detray::actsvg_visualization