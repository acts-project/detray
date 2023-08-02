// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/plugins/actsvg_visualization/svg_converter.hpp"

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

using namespace actsvg;

int main(int, char**) {

    // Generic axes.
    const auto axes =
        draw::x_y_axes("axes", {-250, 250}, {-250, 250}, style::stroke(), "axis1", "axis2");

    // Creating the detector and geomentry context.
    using toy_detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    toy_detector_t::geometry_context context{};

    // Creating the converter for the detector.
    const detray::actsvg_visualization::svg_converter det_converter{det, names};

    // Indexes of the portals in the detector to be visualized.
    std::array indexes{0, 1, 2, 3};

    for (size_t i : indexes) {
        // Getting portal i:
        const auto& description = det.portals()[i];
        const auto portal = detray::surface{det, description};

        // Visualization of portal i:
        const auto svg_xy = det_converter.xy(context, portal);
        detray::actsvg_visualization::write_svg("test_actsvg_portal" + std::to_string(i) + "xy", {axes, svg_xy});
        const auto svg_zr = det_converter.zr(context, portal);
        detray::actsvg_visualization::write_svg("test_actsvg_portal" + std::to_string(i) + "zr", {axes, svg_zr});
    }

}