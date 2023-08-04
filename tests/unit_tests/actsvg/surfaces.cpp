// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/plugins/actsvg_visualization/svg.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/core/defs.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
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

    // Indexes of the surfaces in the detector to be visualized.
    std::array indices{0, 1, 2, 3};

    for (size_t i : indices) {
        std::string name = "test_actsvg_surface" + std::to_string(i);
        // Visualization of portal i:
        const auto svg_xy = converter.xy_surface(name, context, i);
        detray::actsvg_visualization::write_svg(name + "_xy.svg", {axes, svg_xy});
        const auto svg_zr = converter.zr_surface(name, context, i);
        detray::actsvg_visualization::write_svg(name + "_zr.svg", {axes, svg_zr});
    }

}