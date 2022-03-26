/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <matplot/matplot.h>

#include "detray/style/styles.hpp"
#include "detray/tools/generators.hpp"
#include "detray/utils/enumerate.hpp"
#include "detray/view/views.hpp"

namespace detray {

/** Static draw function for arcs
 *
 * @param x0 the x value of the first point
 * @param y0 the y value of the first point
 * @param x1 the x value of the second point
 * @param y1 the y value of the second point
 * @param st is the style class for the grid
 *
 */
static inline void draw_line(scalar x0, scalar y0, scalar x1, scalar y1,
                             const style &st) {
    auto line = matplot::line(x0, y0, x1, y1);
    line->color(st.fill_color);
    line->line_width(st.line_width);
    line->line_style(st.line_style.c_str());
}

/** Static draw function for arcs
 *
 * @param c_x the x value of the center
 * @param c_y the y value of the center
 * @param R the Radius of the arc
 * @param phi_min the min phi value
 * @param phi_max the max phi valu
 * @param st is the style class for the grid
 *
 */
static inline void draw_arc(scalar c_x, scalar c_y, scalar R, scalar phi_min,
                            scalar phi_max, const style &st) {
    matplot::hold(matplot::on);
    std::vector<scalar> theta = matplot::linspace(phi_min, phi_max);
    std::vector<scalar> x = matplot::transform(
        theta, [=](auto theta) { return R * std::cos(theta) + c_x; });
    std::vector<scalar> y = matplot::transform(
        theta, [=](auto theta) { return R * std::sin(theta) + c_y; });
    auto arc = matplot::plot(x, y);
    arc->color(st.fill_color);
    arc->line_width(st.line_width);
    arc->line_style(st.line_style.c_str());
}

/** Static draw function for a polygon
 *
 * @param polygon the polygon points
 * @param st is the style class for the grid
 *
 */
static inline void draw_polygon(const std::vector<point2> &polygon,
                                const style &st) {
    for (size_t j = 1; j <= polygon.size(); ++j) {
        size_t i = j - 1;
        size_t jc = (j == polygon.size()) ? 0 : j;
        auto line = matplot::line(polygon[i][0], polygon[i][1], polygon[jc][0],
                                  polygon[jc][1]);
        line->color(st.fill_color);
        line->line_width(st.line_width);
        line->line_style(st.line_style.c_str());
    }
}

/** Static draw function for vertices
 *
 * @param vertices the vertices generated from the mask
 * @param tf is the transform where the surface is placed
 * @param st is the style class for the grid
 * @param view is the view type for the display
 * @param phi_wrappping check for phi_wrappping in a z-phi view
 *                      (limited to 4 vertices)
 *
 */
template <typename vertex_container_type, typename view_t = single_view>
static inline void draw_vertices(const vertex_container_type &vertices,
                                 const transform3 &tf, const style &st,
                                 const view_t &view = view_t(),
                                 bool phi_wrapping = false) {

    // Create a view of the vertices
    auto [x, y] = view(vertices, tf);

    if (phi_wrapping and vertices.size() == 4) {

        // @note this assumes ordering of the x,y as given by the mask
        if (std::abs(y[0] - y[1]) > M_PI) {

            std::vector<scalar> x_up = {x[0], 0.5 * (x[0] + x[1]),
                                        0.5 * (x[1] + x[2]),
                                        0.5 * (x[2] + x[3]), x[3]};
            std::vector<scalar> y_up = {y[0], M_PI, 2 * M_PI - y[0], M_PI,
                                        y[3]};
            auto filled_area_up = matplot::fill(x_up, y_up, "w");
            filled_area_up->color(st.fill_color);
            filled_area_up->line_width(st.line_width);

            std::vector<scalar> x_down = {x[1], 0.5 * (x[0] + x[1]),
                                          0.5 * (x[1] + x[2]),
                                          0.5 * (x[2] + x[3]), x[2]};
            std::vector<scalar> y_down = {y[1], -2 * M_PI - y[1], -M_PI,
                                          -2 * M_PI - y[2], y[2]};
            auto filled_area_down = matplot::fill(x_down, y_down, "w");
            filled_area_down->color(st.fill_color);
            filled_area_down->line_width(st.line_width);
            return;
        }
    }

    auto filled_area = matplot::fill(x, y, "w");

    filled_area->color(st.fill_color);
    filled_area->line_width(st.line_width);
}

/** Static draw function for masks
 *
 * @param mask is the mask of the surface to be drawn
 * @param tf is the transform where the surface is placed
 * @param st is the style class for the grid
 * @param view is the view type for the display
 */
template <typename mask_type, typename view_t = single_view>
static inline void draw_mask(const mask_type &mask, const transform3 &tf,
                             const style &st, const view_t &view = view_t()) {
    auto contour = vertices(mask, st.segments);
    // Create a view of the vertices
    auto [x, y] = view(contour, tf);
    auto filled_area = matplot::fill(x, y, "w");

    filled_area->color(st.fill_color);
    filled_area->line_width(st.line_width);
}

/** Static draw function for masks
 *
 * @param grid is the grid of the surface to be drawn
 * @param st is the style class for the grid
 * @param view is the view type for the display
 */
template <typename grid_type>
static inline void draw_r_phi_grid(const grid_type &grid, const style &st) {
    const auto &r_axis = grid.axis_p0();
    auto r_borders = r_axis.all_borders();
    for (auto r : r_borders) {
        auto r_line = matplot::rectangle(-r, -r, 2 * r, 2 * r, 1.);
        r_line->color(st.fill_color);
        r_line->line_width(st.line_width);
    }

    const auto &phi_axis = grid.axis_p1();
    auto phi_borders = phi_axis.all_borders();
    for (auto [i, phi] : enumerate(phi_borders)) {
        scalar cos_phi = std::cos(phi);
        scalar sin_phi = std::sin(phi);
        scalar x0 = r_borders[0] * cos_phi;
        scalar y0 = r_borders[0] * sin_phi;
        scalar x1 = r_borders[r_borders.size() - 1] * cos_phi;
        scalar y1 = r_borders[r_borders.size() - 1] * sin_phi;
        if (i < static_cast<size_t>(phi_borders.size() - 1)) {
            auto phi_line = matplot::line(x0, y0, x1, y1);
            phi_line->color(st.fill_color);
            phi_line->line_width(st.line_width);
        }
    }
}

/** Static draw function for masks
 *
 * @param grid is the grid of the surface to be drawn
 * @param st is the style class for the grid
 * @param view is the view type for the display
 */
template <typename grid_type>
static inline void draw_z_phi_grid(const grid_type &grid, const style &st) {
    const auto &z_axis = grid.axis_p0();
    auto z_borders = z_axis.all_borders();
    for (auto z : z_borders) {

        scalar z0 = z;
        scalar phi0 = -M_PI;
        scalar z1 = z;
        scalar phi1 = M_PI;
        auto z_line = matplot::line(z0, phi0, z1, phi1);
        z_line->color(st.fill_color);
        z_line->line_width(st.line_width);
    }

    const auto &phi_axis = grid.axis_p1();
    auto phi_borders = phi_axis.all_borders();
    for (auto [i, phi] : enumerate(phi_borders)) {
        scalar z0 = z_borders[0];
        scalar z1 = z_borders[z_borders.size() - 1];
        auto phi_line = matplot::line(z0, phi, z1, phi);
        phi_line->color(st.fill_color);
        phi_line->line_width(st.line_width);
    }
}

/** Set up a new display with dimensions and potentially wuite
 *
 * @param onscreen The boolean to make it quite or not
 * @param dimensions The dimensions of the supbplot
 *
 * @return an axis object
 */
auto display(bool onscreen = true,
             std::array<float, 4> dimensions = {0.1, 0.1, 0.8, 0.8}) {
    auto ax = matplot::subplot(dimensions);
    ax->parent()->quiet_mode(not onscreen);
    return ax;
}

/** Save the current picture
 *
 * @param filename the name of the file (obviously)
 * @param clearpad an indication if to clear the pad
 */
void save(const std::string &filename, bool clearpad = true) {
    matplot::save(filename);
    if (clearpad) {
        matplot::cla();
    }
}

}  // namespace detray