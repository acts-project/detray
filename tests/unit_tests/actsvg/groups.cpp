// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/actsvg_visualization/svg.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/core/defs.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/meta.hpp"

// System include(s)
#include <array>
#include <string>
#include <vector>

int main(int, char**) {

    // Axes.
    const auto axes =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250}, actsvg::style::stroke(), "axis1", "axis2");

    // Creating the detector and geomentry context.
    using toy_detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    toy_detector_t::geometry_context context{};

    // Creating the converter for the detector.
    const detray::actsvg_visualization::svg_converter converter{det, names};

    // Visualisation of a group of surfaces.
    const std::array surface_group_indices{1, 100, 10, 200};

    const auto svg_surface_group_xy = converter.xy_surfaces("my_surface_group1_xy", context, surface_group_indices);
    detray::actsvg_visualization::write_svg("test_actsvg_surface_group_xy.svg", {axes, svg_surface_group_xy});

    const auto svg_surface_group_zr = converter.zr_surfaces("my_surface_group1_zr", context, surface_group_indices);
    detray::actsvg_visualization::write_svg("test_actsvg_surface_group_zr.svg", {axes, svg_surface_group_zr});
    
    // Visualisation of a group of volumes.
    const std::array volume_group_indices{3, 5};

    const auto svg_volume_group_xy = converter.xy_volumes("my_volume_group1_xy", context, volume_group_indices);
    detray::actsvg_visualization::write_svg("test_actsvg_volume_group_xy.svg", {axes, svg_volume_group_xy});

    const auto svg_volume_group_zr = converter.zr_volumes("my_volume_group1_zr", context, volume_group_indices);
    detray::actsvg_visualization::write_svg("test_actsvg_volume_group_zr.svg", {axes, svg_volume_group_zr});
    
    // Writing SVGs to a combined file.
    // NOTE: The all svg object's identification must be unique in the file!
    detray::actsvg_visualization::write_svg("test_actsvg_volume_and_surface_group.svg", {axes, svg_surface_group_xy, svg_volume_group_zr});
}