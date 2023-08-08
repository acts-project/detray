#pragma once

// Project include(s)
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/proto/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/portal.hpp"

// System include(s)
#include <vector>

namespace detray::actsvg_visualization::proto {

/// @brief Calculates the proto volume of a collection of detray surfaces.
///
/// @param context The context.
/// @param d_surfaces The detray surfaces.
///
/// @returns An actsvg proto volume representing the volume.
template <typename detector_t>
auto volume(
const typename detector_t::geometry_context& context,
const detector_t& detector,
const std::vector<detray::surface<detector_t>>& d_surfaces) {
    proto_volume p_volume;
    std::vector<proto_surface> surfaces;
    std::vector<proto_portal> portals;
    for (const auto& item : d_surfaces){
        if (item.is_portal()){
            auto portal = proto::portal(context, detector, item);
            portals.push_back(portal);
        }
        else{
            auto surface = proto::surface(context, item);
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
template <typename detector_t>
auto volume(
const typename detector_t::geometry_context& context,
const detector_t& detector,
const detray::detector_volume<detector_t>& d_volume) {
    std::vector<detray::surface<detector_t>> surfaces{};
    const auto descriptors = d_volume.surface_lookup();
    for (const auto& desc : descriptors){
        surfaces.push_back(detray::surface{detector, desc});
    }
    return volume(context, detector, surfaces);
}

}