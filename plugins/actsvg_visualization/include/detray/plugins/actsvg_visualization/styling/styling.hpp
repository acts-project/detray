#pragma once

// Project include(s)
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <cstdlib>
#include <vector>

namespace detray::actsvg_visualization::styling {

struct style{
    const std::vector<actsvg::style::color> surface_fill_colors;
    const std::vector<actsvg::style::color> portal_fill_colors;
    const double stroke_width;
    const double marker_size;
    const bool hide_links;
};

namespace colors{

/// @brief Picks a random element in the container.
template <typename container_t>
auto pick_random(container_t container){
    int idx = std::rand() % container.size();
    return container[idx];
}

// Red tones
constexpr std::array cardinal{167,51,63};
constexpr std::array madder{165,28,48};
constexpr std::array auburn{167,51,63};
constexpr std::array burgundy{116,18,29};
constexpr std::array chocolate_cosmos{88,12,31};

// Blue tones
constexpr std::array celestial_blue{62,146,204};
constexpr std::array lapis_lazulli{42,98,143};
constexpr std::array prussian_blue1{19,41,61};
constexpr std::array prussian_blue2{22,50,79};
constexpr std::array indigo_dye{24,67,90};

// Green tones
constexpr std::array celadon{190,230,206};
constexpr std::array aquamarine1{188,255,219};
constexpr std::array aquamarine2{141,255,205};
constexpr std::array emerald{104,216,155};
constexpr std::array shamrock_green{79,157,105};

std::vector<actsvg::style::color> red_theme(const actsvg::scalar opacity){
    return {{cardinal, opacity}, {madder, opacity}, {auburn, opacity}, {burgundy, opacity}, {chocolate_cosmos, opacity}};
}

std::vector<actsvg::style::color> blue_theme(const actsvg::scalar opacity){
    return {{celestial_blue, opacity}, {lapis_lazulli, opacity}, {prussian_blue1, opacity}, {prussian_blue2, opacity}, {indigo_dye, opacity}};
}

std::vector<actsvg::style::color> green_theme(const actsvg::scalar opacity){
    return {{celadon, opacity}, {aquamarine1, opacity}, {aquamarine2, opacity}, {emerald, opacity}, {shamrock_green, opacity}};
}

};

const style default_style{
    colors::blue_theme(0.8),
    colors::red_theme(0.8),
    1.,
    1.2,
    false
};

/// @brief Sets the style of the proto surface.
template <typename point3_container_t>
void apply_style(actsvg::proto::surface<point3_container_t>& p_surface, const style& styling)
{
    auto fill_color = colors::pick_random(styling.surface_fill_colors);
    p_surface._fill = actsvg::style::fill(fill_color);
    p_surface._stroke = actsvg::style::stroke(fill_color, styling.stroke_width);
    
}

/// @brief Sets the style of the proto link.
template <typename point3_container_t>
void apply_style(typename actsvg::proto::portal<point3_container_t>::link& p_link, const style& styling)
{
    p_link._end_marker._size = styling.marker_size;
}

/// @brief Sets the style of the proto portal.
template <typename point3_container_t>
void apply_style(actsvg::proto::portal<point3_container_t>& p_portal, const style& styling)
{
    auto fill_color = colors::pick_random(styling.portal_fill_colors);
    p_portal._surface._fill = actsvg::style::fill(fill_color);
    p_portal._surface._stroke = actsvg::style::stroke(fill_color, styling.stroke_width);
    if (styling.hide_links){
        p_portal._volume_links = {};
    }
    for (auto& volume_link : p_portal._volume_links)
    {
        apply_style<point3_container_t>(volume_link, styling);
    }
}

/// @brief Sets the style of the proto volume.
template <typename point3_container_t>
void apply_style(actsvg::proto::volume<point3_container_t>& p_volume, const style& styling)
{
    for (auto& p_surface : p_volume._v_surfaces){
        apply_style(p_surface, styling);
    }
    for (auto& p_portals : p_volume._portals){
        apply_style(p_portals, styling);
    }
}

} // namespace name