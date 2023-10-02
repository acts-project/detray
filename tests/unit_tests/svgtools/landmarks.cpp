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
#include <string>

int main(int, char**) {

    // Axes.
    const auto axes =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                               actsvg::style::stroke(), "axis1", "axis2");

    // Creating the view.
    const actsvg::views::x_y xy;
    const actsvg::views::z_r zr;

    // Creating the detector and geomentry context.
    using detector_t = detray::detector<detray::toy_metadata>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr);
    detector_t::geometry_context context{};

    using point = typename detector_t::point3;

    // Creating the illustrator class.
    const detray::svgtools::illustrator il{det, context};

    // Sometimes its useful to be able to just draw a point while debugging.
    // For this the draw_landmark function is available.
    const point test_point{100, 50, 20};

    const auto svg_xy = il.draw_landmark("landmark", test_point, xy);
    const auto svg_zr = il.draw_landmark("landmark", test_point, zr);
    detray::svgtools::write_svg("test_svgtools_landmark.svg",
                                {svg_xy, svg_zr, axes});
}