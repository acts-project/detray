#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/proto/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/link.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/link_utils.hpp"
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization::proto {

/// @returns An actsvg proto portal representing the portal.
/// @note detray portal is_portal() must be true.
template <typename detector_t>
auto portal(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_portal)
{
    assert(d_portal.is_portal());
    proto_portal p_portal;
    if (utils::has_link(d_portal))
    {
        p_portal._volume_links = {proto::link(context, d_portal, 3.), proto::link(context, d_portal, -3.)};
    }
    p_portal._surface = surface(context, d_portal);
    return p_portal;
}

}  // namespace detray::actsvg_visualization