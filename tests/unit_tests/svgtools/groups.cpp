/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/plugins/svgtools/utils/groups.hpp"

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
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr);

    // Creating the svg generator for the detector.
    const detray::svgtools::illustrator il{det, names};

    // Visualisation of a group of surfaces.
    const std::array surface_group_indices{1UL, 100UL, 10UL, 200UL};

    const auto svg_surface_group_xy =
        il.draw_surfaces("my_surface_group1_xy", surface_group_indices, xy);
    detray::svgtools::write_svg("test_svgtools_surface_group_xy",
                                {axes, svg_surface_group_xy});

    const auto svg_surface_group_zr =
        il.draw_surfaces("my_surface_group1_zr", surface_group_indices, zr);
    detray::svgtools::write_svg("test_svgtools_surface_group_zr.svg",
                                {axes, svg_surface_group_zr});

    // Visualisation of a group of volumes.
    const std::array volume_group_indices{3UL, 5UL};

    const auto svg_volume_group_xy =
        il.draw_volumes("my_volume_group1_xy", volume_group_indices, xy);
    detray::svgtools::write_svg("test_svgtools_volume_group_xy",
                                {axes, svg_volume_group_xy});

    const auto svg_volume_group_zr =
        il.draw_volumes("my_volume_group1_zr", volume_group_indices, zr);
    detray::svgtools::write_svg("test_svgtools_volume_group_zr",
                                {axes, svg_volume_group_zr});

    // We can also use the svgtools::utils to group actsvg::svg::objects into
    // one.
    auto svg_combined_group = detray::svgtools::utils::group("combined_group");
    svg_combined_group.add_object(svg_surface_group_xy);
    svg_combined_group.add_object(svg_volume_group_zr);

    detray::svgtools::write_svg("test_svgtools_combined_group1.svg",
                                {axes, svg_combined_group});

    // Alternatively, this is equivalent to:
    detray::svgtools::write_svg(
        "test_svgtools_combined_group2.svg",
        {axes, svg_surface_group_xy, svg_volume_group_zr});

    // NOTE: The all svg object's identification must be unique in the
    // file!
}
