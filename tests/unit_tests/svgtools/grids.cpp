/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

int main(int, char**) {

    // This test creates the visualization using the illustrator class.
    // However, for full control over the process, it is also possible to use
    // the tools in svgstools::conversion, svgstools::display, and
    // actsvg::display by converting the object to a proto object, optionally
    // styling it, and then displaying it.

    // Creating the detector and geomentry context.
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr);
    using detector_t = decltype(det);
    detector_t::geometry_context context{};

    // Creating the view.
    const actsvg::views::x_y view;

    // Creating the svg generator for the detector.
    detray::svgtools::illustrator il{det, context, true};

    // In this example we want to draw the grids of the volumes with indices 0,
    // 1, ... in the detector.
    std::vector<std::size_t> indices = {
        0UL,  1UL,  2UL,  3UL,  4UL,  5UL,  6UL,  7UL,  8UL,  9UL,
        10UL, 11UL, 12UL, 13UL, 14UL, 15UL, 16UL, 17UL, 18UL, 19UL};

    for (const auto i : indices) {
        std::string name = "volume" + std::to_string(i) + "_grid";
        // Draw the grid for volume i.
        if (auto grid_svg_ptr = il.draw_grid(name, i, view)) {
            // Draw volume i.
            const auto volume_svg = il.draw_volume("volume", i, view);
            // Write volume i and its grid
            // (the grid is last in the list since we want it to be displayed on
            // top of the volume).
            detray::svgtools::write_svg(name + ".svg",
                                        {volume_svg, *grid_svg_ptr});
        }
    }
}