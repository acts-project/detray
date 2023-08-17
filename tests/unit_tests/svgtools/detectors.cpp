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

    actsvg::style::stroke stroke_black = actsvg::style::stroke();
    // x-y axis.
    auto xy_axis =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250}, stroke_black, "x", "y");  
    // z-r axis.
    auto zr_axis =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250}, stroke_black, "z", "r");  

    // Creating the views.
    const actsvg::views::x_y xy;
    const actsvg::views::z_r zr;

    // Creating the detector and geomentry context.
    using toy_detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    toy_detector_t::geometry_context context{};

    // Creating the svg generator for the detector.
    detray::svgtools::illustrator il{det, names};

    // Get the svgtools of the toy detetector in x-y view.
    const auto svg_xy = il.draw_detector("detector_xy", context, xy);
    // Write the svgtools of toy detector.
    detray::svgtools::write_svg("test_svgtools_detector_xy.svgtools", {xy_axis, svg_xy});

    // Get the svgtools of the toy detetector in z-r view.
    const auto svg_zr = il.draw_detector("detector_zr", context, zr);
    // Write the svgtools of toy detector in z-r view
    detray::svgtools::write_svg("test_svgtools_detector_zr.svgtools", {zr_axis, svg_zr});

}