#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/plugins/actsvg/portal_conversion.hpp"
#include "detray/plugins/actsvg/surface_conversion.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/volume.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

using point3 = std::array<actsvg::scalar, 3>;
using point3_container = std::vector<point3>;
using proto_surface = actsvg::proto::surface<point3_container>;
using proto_portal = actsvg::proto::portal<point3_container>;
using proto_volume = actsvg::proto::volume<point3_container>;

/// @brief Calculates the proto volume of a detray volume.
///
/// @param d_volume The detray volume.
/// @param context The context.
///
/// @returns An actsvg proto volume representing the volume.
template <typename detector_t>
auto convert_volume(
    const detector_t& detector,
    const detray::detector_volume<detector_t>& d_volume,
    const typename detector_t::geometry_context& context) {
        proto_volume p_volume;
        std::vector<proto_surface> surfaces;
        std::vector<proto_portal> portals;
        for (auto& description : d_volume.surface_lookup()){
            const auto item = detray::surface{detector, description};
            if (item.is_portal()){
                auto portal = convert_portal(detector, item, context);
                portals.push_back(portal);
            }
            else{
                auto surface = convert_surface(item, context);
                surfaces.push_back(surface);
            }
        }
        p_volume._v_surfaces = surfaces;
        p_volume._portals = portals;
        return p_volume;
}

}