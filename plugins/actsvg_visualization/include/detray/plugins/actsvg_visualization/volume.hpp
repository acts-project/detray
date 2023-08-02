#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/plugins/actsvg_visualization/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/surface.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/volume.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization::volume {

struct volume_options{
    surface::surface_options s_options;
    surface::portal_options p_options;
    // Indexes of the visible surfaces/portals.
    std::vector<int> visible_surfaces;
};

/// @brief Calculates the proto volume of a detray volume.
///
/// @param d_volume The detray volume.
/// @param context The context.
///
/// @returns An actsvg proto volume representing the volume.
template <typename detector_t>
auto to_proto_volume(
const typename detector_t::geometry_context& context,
const detector_t& detector,
const detray::detector_volume<detector_t>& d_volume,
const volume_options& v_options) {
    conversion_types::volume p_volume;
    std::vector<conversion_types::surface> surfaces;
    std::vector<conversion_types::portal> portals;
    for (auto& description : d_volume.surface_lookup()){
        const auto item = detray::surface{detector, description};
        if (item.is_portal()){
            auto portal = surface::to_proto_portal(context, item, v_options.p_options);
            portals.push_back(portal);
        }
        else{
            auto surface = surface::to_proto_surface(context, item, v_options.s_options);
            surfaces.push_back(surface);
        }
    }
    p_volume._v_surfaces = surfaces;
    p_volume._portals = portals;
    return p_volume;
}

template <typename detector_t, typename view_t>
auto to_svg(const typename detector_t::geometry_context& context, const view_t& view, const detector_t& detector, const detray::detector_volume<detector_t>& d_volume, const volume_options& v_options, const std::string& name)
{
    auto p_volume = volume::to_proto_volume(context, detector, d_volume, v_options);
    return actsvg::display::volume(name, p_volume, view);
}
}