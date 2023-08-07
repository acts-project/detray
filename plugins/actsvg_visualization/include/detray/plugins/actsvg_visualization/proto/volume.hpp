#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/proto/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/portal.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/volume.hpp"

// System include(s)
#include <type_traits>
#include <vector>

#include <iostream>

namespace detray::actsvg_visualization::proto {

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
    proto_volume p_volume;
    std::vector<proto_surface> surfaces;
    std::vector<proto_portal> portals;
    for (const auto desc : d_volume.surface_lookup()){
        const auto item = detray::surface{detector, desc};
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
            auto surface = proto::surface(context, detector, item);
            surfaces.push_back(surface);
        }
    }
    p_volume._v_surfaces = surfaces;
    p_volume._portals = portals;
    return p_volume;
}

}