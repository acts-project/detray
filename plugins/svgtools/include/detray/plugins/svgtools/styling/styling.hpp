/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/plugins/svgtools/meta/proto/eta_lines.hpp"
#include "detray/plugins/svgtools/meta/proto/intersection.hpp"
#include "detray/plugins/svgtools/meta/proto/landmark.hpp"
#include "detray/plugins/svgtools/meta/proto/trajectory.hpp"
#include "detray/plugins/svgtools/styling/colors.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/core/style.hpp"
#include "actsvg/proto/detector.hpp"
#include "actsvg/proto/grid.hpp"
#include "actsvg/proto/portal.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"
#include "actsvg/styles/defaults.hpp"

// System include(s)
#include <cstdlib>
#include <vector>

namespace detray::svgtools::styling {

/// Style applied to an actsvg proto surface
struct surface_style {

    surface_style(actsvg::style::color fill_color =
                      colors::tableau_colorblind10::red_tones(0.4f).front(),
                  actsvg::scalar stroke_width = 0.4f)
        : _fill_color{fill_color}, _stroke_width{stroke_width} {

        _stroke_color = fill_color;
        _highlight_rgb = colors::macaroni_and_cheese;
        _highlights = {};
        _highlight_stroke_rgb = _highlight_rgb;
        _highlight_stroke_width = stroke_width;
        _font_size = 14u;
        _n_segments = 72u;
    }
    // Fill color
    actsvg::style::color _fill_color;

    // Stroke
    actsvg::scalar _stroke_width;
    actsvg::style::color _stroke_color;

    // Highlights
    std::array<int, 3> _highlight_rgb;
    std::vector<std::string> _highlights;
    std::array<int, 3> _highlight_stroke_rgb;
    actsvg::scalar _highlight_stroke_width;

    unsigned int _font_size;
    // Granularity of vertex approximation of arcs
    unsigned int _n_segments;
};

/// Style applied to an actsvg proto portal link
struct link_style {
    actsvg::scalar _marker_size;
};

/// Style applied to an actsvg proto portal
struct portal_style {
    surface_style _surface_style;
    link_style _link_style;
};

/// Style applied to an actsvg proto grid
struct grid_style {
    actsvg::style::color _stroke_color;
    actsvg::scalar _stroke_width;
};

/// Style applied to an actsvg proto volume
struct volume_style {
    surface_style _surface_style;
    portal_style _portal_style;
    grid_style _grid_style;
};

/// Style applied to an actsvg proto detector
struct detector_style {
    volume_style _volume_style;
    grid_style _grid_style;
};

/// Style applied to proto eta lines
struct eta_lines_style {
    actsvg::style::color _fill_color_main;
    actsvg::style::color _fill_color_half;
    actsvg::scalar _stroke_width_main;
    actsvg::scalar _stroke_width_half;
};

/// Style applied to a proto landmark
struct landmark_style {
    actsvg::style::color _fill_color;
    actsvg::scalar _stroke_width;
    actsvg::scalar _marker_size;
    std::string _marker_type;
};

/// Style applied to a proto trajectory
struct trajectory_style {
    // Give different tracks different colors
    actsvg::style::color _fill_color;
    actsvg::scalar _stroke_width;
};

/// Global styling options
struct style {
    detector_style _detector_style;
    eta_lines_style _eta_lines_style;
    trajectory_style _trajectory_style;
    landmark_style _landmark_style;
    landmark_style _intersection_style;
};

/// The default styling
namespace svg_default {

// Surface styles
const styling::surface_style surface_style{colors::blue_theme(0.5f).front(),
                                           3.f};
const styling::surface_style surface_style_portal{
    colors::red_theme(0.5f).front(), 3.f};

// Portal style
const styling::link_style link_style{1.2f};
const styling::portal_style portal_style{surface_style_portal, link_style};

// Volume style
const styling::grid_style grid_style{colors::red_theme(1.f).front(), 1.f};
const styling::volume_style volume_style{surface_style, portal_style,
                                         grid_style};

// Detector style
const styling::detector_style detector_style{volume_style, grid_style};

// Eta lines style
const styling::eta_lines_style eta_lines_style{
    {colors::black, 1.f}, {colors::black, 1.f}, 1.f, 1.f};

// Intersection points and track state style
const styling::landmark_style intersection_style{
    colors::black_theme(1.f).front(), 0.8f, 5.f, "x"};
const styling::landmark_style landmark_style{colors::green_theme(1.f).front(),
                                             0.8f, 3.f, "o"};

// Trajectory style
const styling::trajectory_style trajectory_style{
    colors::green_theme(1.f).front(), 1.f};

// Full style
const styling::style style{detector_style, eta_lines_style, trajectory_style,
                           landmark_style, intersection_style};
}  // namespace svg_default

/// Styling that matches data plotting
namespace tableau_colorblind {

// Surface styles
const styling::surface_style surface_style_passive{
    colors::tableau_colorblind10::grey_tones(0.8f).front(), 1.f};
const styling::surface_style surface_style_portal{
    colors::tableau_colorblind10::blue_tones(0.2f).front(), 1.5f};
const styling::surface_style surface_style_sensitive{
    colors::tableau_colorblind10::red_tones(0.3f).front(), 1.f};

// Portal style
const styling::portal_style portal_style{surface_style_portal,
                                         svg_default::link_style};

// Volume style
const styling::grid_style grid_style{{colors::red, 1.f}, 0.6f};
const styling::volume_style volume_style{surface_style_sensitive, portal_style,
                                         grid_style};

// Detector style
const styling::detector_style detector_style{volume_style, grid_style};

// Full style
const styling::style style{
    detector_style, svg_default::eta_lines_style, svg_default::trajectory_style,
    svg_default::landmark_style, svg_default::intersection_style};
}  // namespace tableau_colorblind

/// @brief Sets the style of the proto surface.
template <typename point3_container_t>
inline void apply_style(actsvg::proto::surface<point3_container_t>& p_surface,
                        const surface_style& styling) {
    // Fill color
    p_surface._fill = actsvg::style::fill(styling._fill_color);
    p_surface._fill._fc._hl_rgb = styling._highlight_rgb;
    p_surface._fill._fc._highlight = styling._highlights;

    // Stroke
    p_surface._stroke =
        actsvg::style::stroke(styling._stroke_color, styling._stroke_width);
    p_surface._stroke._sc._hl_rgb = styling._highlight_stroke_rgb;
    p_surface._stroke._hl_width = styling._highlight_stroke_width;
}

/// @brief Sets the style of the proto link.
template <typename point3_container_t>
inline void apply_style(
    typename actsvg::proto::portal<point3_container_t>::link& p_link,
    const link_style& styling) {
    p_link._end_marker._size = styling._marker_size;
}

/// @brief Sets the style of the proto portal.
template <typename point3_container_t>
inline void apply_style(actsvg::proto::portal<point3_container_t>& p_portal,
                        const portal_style& styling) {

    apply_style(p_portal._surface, styling._surface_style);

    for (auto& volume_link : p_portal._volume_links) {
        apply_style<point3_container_t>(volume_link, styling._link_style);
    }
}

/// @brief Sets the style of the proto grid.
inline void apply_style(actsvg::proto::grid& p_grid,
                        const grid_style& styling) {
    p_grid._stroke._sc = styling._stroke_color;
    p_grid._stroke._width = styling._stroke_width;
    // Use dashed lines for the grid
    // p_grid._stroke._dasharray = {4};
}

/// @brief Sets the style of the proto volume.
template <typename point3_container_t>
inline void apply_style(actsvg::proto::volume<point3_container_t>& p_volume,
                        const volume_style& styling) {
    for (auto& p_surface : p_volume._v_surfaces) {
        apply_style(p_surface, styling._surface_style);
    }
    for (auto& p_portals : p_volume._portals) {
        apply_style(p_portals, styling._portal_style);
    }
    apply_style(p_volume._surface_grid, styling._grid_style);
}

/// @brief Sets the style of the proto volume.
template <typename point3_container_t>
inline void apply_style(actsvg::proto::detector<point3_container_t>& p_detector,
                        const detector_style& styling) {
    for (auto& p_volume : p_detector._volumes) {
        apply_style(p_volume, styling._volume_style);
    }
}

/// @brief Sets the style of the proto eta_lines.
inline void apply_style(meta::proto::eta_lines& p_eta_lines,
                        const eta_lines_style& styling) {
    p_eta_lines._stroke_main = actsvg::style::stroke(
        styling._fill_color_main, styling._stroke_width_main);

    p_eta_lines._stroke_half = actsvg::style::stroke(
        styling._fill_color_half, styling._stroke_width_half);
    p_eta_lines._stroke_half._dasharray = {2, 2};
}

/// @brief Sets the style of the proto landmark.
template <typename point3_t>
inline void apply_style(meta::proto::landmark<point3_t>& p_landmark,
                        const landmark_style& styling) {

    const auto fill = actsvg::style::fill(styling._fill_color);
    const auto stroke =
        actsvg::style::stroke(styling._fill_color, styling._stroke_width);

    const auto marker = actsvg::style::marker{
        styling._marker_type, styling._marker_size, fill, stroke};
    p_landmark._marker = marker;
}

/// @brief Sets the style of the proto intersection record.
template <typename point3_t>
inline void apply_style(meta::proto::intersection<point3_t>& p_intersection,
                        const landmark_style& styling) {

    const auto stroke =
        actsvg::style::stroke(styling._fill_color, styling._stroke_width);

    auto fill_color = styling._fill_color;
    fill_color._opacity *= 0.5f;
    const auto fill = actsvg::style::fill(fill_color);

    const auto marker = actsvg::style::marker{
        styling._marker_type, styling._marker_size, fill, stroke};

    for (auto& p_landmark : p_intersection._landmarks) {
        p_landmark._marker = marker;
    }
}

/// @brief Sets the style of the proto trajectory.
template <typename point3_t>
inline void apply_style(meta::proto::trajectory<point3_t>& p_trajectory,
                        const trajectory_style& styling) {
    p_trajectory._stroke =
        actsvg::style::stroke(styling._fill_color, styling._stroke_width);
}

}  // namespace detray::svgtools::styling
