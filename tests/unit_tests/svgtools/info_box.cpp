/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/masks/masks.hpp"

#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/writer.hpp"
#include "detray/plugins/svgtools/conversion/information_section.hpp"
#include "detray/plugins/svgtools/meta/display/geometry.hpp"
#include "detray/plugins/svgtools/meta/display/information.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"

// System include(s)
#include <array>
#include <string>
#include <iostream>


int main(int, char**) {

    // Axes.
    const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                             actsvg::style::stroke());

    const actsvg::views::x_y view;

    using point3 = std::array<actsvg::scalar, 3>;
    using point3_container = std::vector<point3>;

    using toy_detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    toy_detector_t::geometry_context context{};

    const auto surface = detray::surface{det, det.surface_lookup()[0]};
    const auto proto_sur = detray::svgtools::conversion::surface<point3_container>(context, surface);
    const auto svg_sur = actsvg::display::surface("surface", proto_sur, view);

    const auto proto_is = detray::svgtools::conversion::information_section<point3>(context, surface);
    const actsvg::point2 offset{100, 100};
    const auto svg_is = detray::svgtools::meta::display::information_section("info", proto_is, view, offset, svg_sur);

    detray::svgtools::write_svg("test_info.svg", {axes, svg_sur, svg_is});

    //TODO: Test on a surface.
    //TODO: Rename or put text code in another file than "geometry" .hpp.
}