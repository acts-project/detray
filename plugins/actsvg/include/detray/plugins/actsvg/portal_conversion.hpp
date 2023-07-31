#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/link_conversion.hpp"
#include "detray/plugins/actsvg/surface_conversion.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/core.hpp"

// Algebra include(s)
#include "algebra/math/cmath.hpp"

// System include(s)
#include <assert.h>
#include <vector>
#include <type_traits>
#include <iostream>

namespace detray::actsvg_visualization {
/// @returns An actsvg proto portal representing the portal.
/// @note detray portal is_portal() must be true.
template <typename detector_t>
proto_portal convert_portal(const detector_t& detector, const detray::surface<detector_t>& d_portal, const typename detector_t::geometry_context& context)
{
    assert(d_portal.is_portal());
    proto_portal p_portal;
    if (has_link(d_portal))
    {
        const auto d_volume = get_link_volume(detector, d_portal);
        const auto p_link = convert_link(d_portal, context);
        p_portal._volume_links = {p_link};
    }
    p_portal._surface = convert_surface(d_portal, context);
    return p_portal;
}
}