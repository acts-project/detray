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

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <array>
#include <string>

int main(int, char**) {

    // Axes.
    const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                             actsvg::style::stroke());

    // Creating the views.
    const actsvg::views::x_y xy;
    const actsvg::views::z_r zr;

    // Creating the detector and geomentry context.
    using toy_detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    toy_detector_t::geometry_context context{};

    // Creating the svg generator for the detector.
    const detray::svgtools::illustrator il{det, names};

    // Indexes of the volumes in the detector to be visualized.
    std::array indices{0, 1, 2, 3};
    for (int i : indices) {
        std::string name = "test_svgtools_volume" + std::to_string(i);
        // Visualization of volume i:
        const auto svg_xy = il.draw_volume(name, context, i, xy);
        detray::svgtools::write_svg(name + "_xy.svg", {axes, svg_xy});
        const auto svg_zr = il.draw_volume(name, context, i, zr);
        detray::svgtools::write_svg(name + "_zr.svg", {axes, svg_zr});
    }
}