/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "detray/plugins/svgtools/meta/proto/intersection_record.hpp"
#include "detray/plugins/svgtools/meta/proto/landmark.hpp"

// System include(s)
#include <cstdlib>
#include <vector>

namespace detray::svgtools::styling {

namespace colors {

/// @brief Picks a random element in the container.
template <typename container_t>
auto pick_random(container_t container) {
    const auto idx = static_cast<std::size_t>(std::rand()) % container.size();
    return container[idx];
}

// Black
constexpr std::array black{0, 0, 0};

// Red tones
constexpr std::array cardinal{167, 51, 63};
constexpr std::array madder{165, 28, 48};
constexpr std::array auburn{167, 51, 63};
constexpr std::array burgundy{116, 18, 29};
constexpr std::array chocolate_cosmos{88, 12, 31};

// Blue tones
constexpr std::array celestial_blue{62, 146, 204};
constexpr std::array lapis_lazulli{42, 98, 143};
constexpr std::array prussian_blue1{19, 41, 61};
constexpr std::array prussian_blue2{22, 50, 79};
constexpr std::array indigo_dye{24, 67, 90};

// Green tones
constexpr std::array celadon{190, 230, 206};
constexpr std::array aquamarine1{188, 255, 219};
constexpr std::array aquamarine2{141, 255, 205};
constexpr std::array emerald{104, 216, 155};
constexpr std::array shamrock_green{79, 157, 105};

std::vector<actsvg::style::color> black_theme(const actsvg::scalar opacity) {
    return {{black, opacity}};
}

std::vector<actsvg::style::color> red_theme(const actsvg::scalar opacity) {
    return {{cardinal, opacity},
            {madder, opacity},
            {auburn, opacity},
            {burgundy, opacity},
            {chocolate_cosmos, opacity}};
}

std::vector<actsvg::style::color> blue_theme(const actsvg::scalar opacity) {
    return {{celestial_blue, opacity},
            {lapis_lazulli, opacity},
            {prussian_blue1, opacity},
            {prussian_blue2, opacity},
            {indigo_dye, opacity}};
}

std::vector<actsvg::style::color> green_theme(const actsvg::scalar opacity) {
    return {{celadon, opacity},
            {aquamarine1, opacity},
            {aquamarine2, opacity},
            {emerald, opacity},
            {shamrock_green, opacity}};
}

}  // namespace colors

struct grid_style {
    actsvg::scalar _stroke_width;
};

struct landmark_style {
    std::vector<actsvg::style::color> _fill_colors;
    actsvg::scalar _marker_size;
    std::string _marker_type;
};

struct trajectory_style {
    std::vector<actsvg::style::color> _fill_colors;
    actsvg::scalar _stroke_width;
};

struct surface_style {
    std::vector<actsvg::style::color> _fill_colors;
    actsvg::scalar _stroke_width;
};

struct link_style {
    actsvg::scalar _marker_size;
};

struct portal_style {
    surface_style _surface_style;
    link_style _link_style;
    bool _hide_links;
};

struct volume_style {
    surface_style _surface_style;
    portal_style _portal_style;
};

struct style {
    volume_style _volume_style;
    landmark_style _intersection_style;
    trajectory_style _trajectory_style;
    grid_style _grid_style;
    landmark_style _landmark_style;
};

const surface_style surface_style1{colors::blue_theme(0.5f), 3.f};

const surface_style surface_style2{colors::red_theme(0.5f), 3.f};

const link_style link_style1{1.2f};

const portal_style portal_style1{surface_style2, link_style1, false};

const volume_style volume_style1{surface_style1, portal_style1};

const landmark_style landmark_style1{colors::black_theme(0.8f), 5.f, "x"};

const landmark_style landmark_style2{colors::black_theme(1.f), 3.f, "o"};

const grid_style grid_style1{2.f};

const trajectory_style trajectory_style1{colors::green_theme(1.f), 3.f};

const style style1{volume_style1, landmark_style1, trajectory_style1, grid_style1, landmark_style2};

/// @brief Sets the style of the proto surface.
template <typename point3_container_t>
void apply_style(actsvg::proto::surface<point3_container_t>& p_surface,
                 const surface_style& styling) {
    auto fill_color = colors::pick_random(styling._fill_colors);
    p_surface._fill = actsvg::style::fill(fill_color);
    p_surface._stroke =
        actsvg::style::stroke(fill_color, styling._stroke_width);
}

/// @brief Sets the style of the proto link.
template <typename point3_container_t>
void apply_style(
    typename actsvg::proto::portal<point3_container_t>::link& p_link,
    const link_style& styling) {
    p_link._end_marker._size = styling._marker_size;
}

/// @brief Sets the style of the proto portal.
template <typename point3_container_t>
void apply_style(actsvg::proto::portal<point3_container_t>& p_portal,
                 const portal_style& styling) {
    auto fill_color = colors::pick_random(styling._surface_style._fill_colors);
    p_portal._surface._fill = actsvg::style::fill(fill_color);
    p_portal._surface._stroke =
        actsvg::style::stroke(fill_color, styling._surface_style._stroke_width);
    if (styling._hide_links) {
        p_portal._volume_links = {};
    }
    for (auto& volume_link : p_portal._volume_links) {
        apply_style<point3_container_t>(volume_link, styling._link_style);
    }
}

/// @brief Sets the style of the proto volume.
template <typename point3_container_t>
void apply_style(actsvg::proto::volume<point3_container_t>& p_volume,
                 const volume_style& styling) {
    for (auto& p_surface : p_volume._v_surfaces) {
        apply_style(p_surface, styling._surface_style);
    }
    for (auto& p_portals : p_volume._portals) {
        apply_style(p_portals, styling._portal_style);
    }
}

/// @brief Sets the style of the proto landmark.
template <typename point3_t>
void apply_style(
    meta::proto::landmark<point3_t>& p_landmark,
    const landmark_style& styling) {
    const auto fill_color = colors::pick_random(styling._fill_colors);
    const auto fill = actsvg::style::fill(fill_color);
    const auto stroke = actsvg::style::stroke(fill_color, 3.f);
    const auto marker =
        actsvg::style::marker{styling._marker_type, styling._marker_size, fill, stroke};
    p_landmark._marker = marker;
}


/// @brief Sets the style of the proto intersection record.
template <typename point3_t>
void apply_style(
    meta::proto::intersection_record<point3_t>& p_intersection_record,
    const landmark_style& styling) {
    const auto fill_color = colors::pick_random(styling._fill_colors);
    const auto fill = actsvg::style::fill(fill_color);
    const auto stroke = actsvg::style::stroke(fill_color, 3.f);
    const auto marker =
        actsvg::style::marker{styling._marker_type, styling._marker_size, fill, stroke};
    for (auto& p_landmark : p_intersection_record._landmarks) {
        p_landmark._marker = marker;
    }
}

/// @brief Sets the style of the proto trajectory.
template <typename point3_t>
void apply_style(meta::proto::trajectory<point3_t>& p_trajectory,
                 const trajectory_style& styling) {
    auto fill_color = colors::pick_random(styling._fill_colors);
    p_trajectory._stroke =
        actsvg::style::stroke(fill_color, styling._stroke_width);
}

/// @brief Sets the style of the proto grid.
void apply_style(actsvg::proto::grid& p_grid,
                 const grid_style& styling) {
    p_grid._stroke._width = styling._stroke_width;
}

template <typename style1_t, typename style2_t>
auto copy_fill_colors(style1_t target, const style2_t& reference) {
    target._fill_colors = reference._fill_colors;
    return target;
}

}  // namespace detray::svgtools::styling
