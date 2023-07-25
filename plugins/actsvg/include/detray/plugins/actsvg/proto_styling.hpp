#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/portal_conversion.hpp"
#include "detray/plugins/actsvg/surface_conversion.hpp"
#include "detray/plugins/actsvg/volume_conversion.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

struct surface_style
{
    const std::vector<actsvg::style::color> fill_colors;
};

struct detector_style
{
    surface_style surface_styling;
    surface_style portal_styling;
};

const surface_style default_surface_style{{
    {{0, 128, 0}, 0.5},
    {{34, 139, 34}, 0.5},
    {{60, 179, 113}, 0.5},
    {{85, 107, 47}, 0.5},
    {{124, 252, 0}, 0.5},
    {{154, 205, 50}, 0.5}
    
    }};
const surface_style default_portal_style{{
    
    {{128, 0, 128}, 0.5},
    {{138, 43, 226}, 0.5},
    {{148, 0, 211}, 0.5},
    {{160, 32, 240}, 0.5},
    {{186, 85, 211}, 0.5},
    {{218, 112, 214}, 0.5}
    
    }};


const detector_style default_detector_style{default_surface_style, default_portal_style};
/// @brief Picks a random element in the container.
template <typename container_t>
auto pick_random(container_t container){
    int idx = rand() % container.size();
    return container[idx];
}

using proto_surface = actsvg::proto::surface<point3_container>;
using proto_portal = actsvg::proto::portal<point3_container>;
using proto_volume = actsvg::proto::volume<point3_container>;

/// @brief Sets the style of the proto surface.
template <typename point3_container_t>
void apply_style(actsvg::proto::surface<point3_container_t>& p_surface, const detector_style& style)
{
    auto fill_color = pick_random(style.surface_styling.fill_colors);
    p_surface._fill = actsvg::style::fill(fill_color);
}

/// @brief Sets the style of the proto portal.
template <typename point3_container_t>
void apply_style(actsvg::proto::portal<point3_container_t>& p_portal, const detector_style& style)
{
    auto fill_color = pick_random(style.portal_styling.fill_colors);
    p_portal._surface._fill = actsvg::style::fill(fill_color);
}

/// @brief Sets the style of the proto volume.
template <typename point3_container_t>
void apply_style(actsvg::proto::volume<point3_container_t>& p_volume, const detector_style& style)
{
    for (auto& p_surface : p_volume._v_surfaces){
        apply_style(p_surface, style);
    }
    for (auto& p_portals : p_volume._portals){
        apply_style(p_portals, style);
    }
}

/// @brief Sets the style of the proto detector.
template <typename point3_container_t>
void apply_style(actsvg::proto::detector<point3_container_t>& p_detector, const detector_style& style)
{
    for (auto& p_volume : p_detector._volumes){
        apply_style(p_volume, style);
    }
}
}  // namespace detray::actsvg_visualization