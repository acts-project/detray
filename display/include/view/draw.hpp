/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "view/views.hpp"
#include "view/generators.hpp"
#include "style/styles.hpp"
#include "utils/enumerate.hpp"

#include <matplot/matplot.h>

namespace detray
{

    /** Static draw function for vertices 
     * 
     * @param vertices the vertices generated from the mask
     * @param tf is the transform where the surface is placed
     * @param st is the style class for the grid
     * @param view is the view type for the display
     * 
     */
    template <typename vertex_container_type, typename view_t = single_view>
    static inline void draw_vertices(const vertex_container_type &vertices,
                                     const transform3 &tf,
                                     const style &st,
                                     const view_t &view = view_t())
    {
        // Create a view of the vertices
        auto [x, y] = view(vertices, tf);
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
    static inline void draw_mask(const mask_type &mask,
                                 const transform3 &tf,
                                 const style &st,
                                 const view_t &view = view_t())
    {
        // Create a view of the vertices
        auto [x, y] = view(vertices(mask, st.segments), tf);
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
    static inline void draw_r_phi_grid(const grid_type &grid,
                                       const style &st)
    {
        const auto &r_axis = grid.axis_p0();
        auto r_borders = r_axis.all_borders();
        for (auto r : r_borders)
        {
            auto r_line = matplot::rectangle(-r, -r, 2 * r, 2 * r, 1.);
            r_line->color(st.fill_color);
            r_line->line_width(st.line_width);
        }

        const auto &phi_axis = grid.axis_p1();
        auto phi_borders = phi_axis.all_borders();
        for (auto [i, phi] : enumerate(phi_borders))
        {
            scalar cos_phi = std::cos(phi);
            scalar sin_phi = std::sin(phi);
            scalar x0 = r_borders[0] * cos_phi;
            scalar y0 = r_borders[0] * sin_phi;
            scalar x1 = r_borders[r_borders.size() - 1] * cos_phi;
            scalar y1 = r_borders[r_borders.size() - 1] * sin_phi;
            if (i < static_cast<size_t>(phi_borders.size() - 1))
            {
                auto phi_bin = matplot::line(x0, y0, x1, y1);
                phi_bin->color(st.fill_color);
                phi_bin->line_width(st.line_width);
            }
        }
    }

    /** Set up a new display with dimensions and potentially wuite
     *
     * @param onscreen The boolean to make it quite or not
     * @param dimensions The dimensions of the supbplot
     * 
     * @return an axis object
     */
    auto display(bool onscreen = true, std::array<float, 4> dimensions = {0.1, 0.1, 0.8, 0.8})
    {
        auto ax = matplot::subplot(dimensions);
        ax->parent()->quiet_mode(not onscreen);
        return ax;
    }

    /** Save the current picture
     *
     * @param filename the name of the file (obviously)
     * @param clearpad an indication if to clear the pad
     */
    void save(const std::string &filename, bool clearpad = true)
    {
        matplot::save(filename);
        if (clearpad)
        {
            matplot::cla();
        }
    }

} // namespace detray