#pragma once

#include "view/views.hpp"
#include "view/generators.hpp"
#include "style/styles.hpp"

#include <matplot/matplot.h>
#include <iostream>

namespace detray
{

    /** Static draw function for the surface (components)
     *   
     * @param maks is the mask of the surface to be drawn
     * @param tf is the transform where the surface is placed
     * @param st is the style class for the grid
     * @param view is the view type for the display
     */
    template <typename mask_type, typename view_t = single_view>
    static inline void draw(const mask_type &mask,
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

    /** Show what's currently to show */
    void show()
    {
        matplot::show();
    }

} // namespace aplot