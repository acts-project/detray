// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/masks/cylinder2D.hpp"
#include "detray/masks/masks.hpp"
#include "detray/plugins/actsvg/svg_conversion.hpp"
#include "detray/tracks/tracks.hpp"

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

    using toy_detector_t = detray::detector<detray::toy_metadata<>>;

    vecmem::host_memory_resource host_mr;
    const toy_detector_t det = detray::create_toy_geometry(host_mr, 4, 3);
    toy_detector_t::geometry_context context{};

    views::x_y view;

    // Create SVG file.
    svg::file file;

    auto indices = {13, 20, 100, 150, 200, 250};

    std::vector<actsvg::svg::object> objs;

    for (const auto& pair :
         detray::views::pick(det.surface_lookup(), indices)) {
        const auto index = pair.first;
        const auto description = pair.second;
        const auto surface = detray::surface{det, description};

        auto name = "detector_surface" + std::to_string(index);
        const auto svg = detray::actsvg_visualization::svg(name, det, surface, context, view);

        objs.push_back(svg);
    }

    // Draw x-y-axis.
    style::stroke stroke_black = style::stroke();
    auto axis =
        draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");  
    objs.push_back(axis);

    detray::actsvg_visualization::write_svg(objs, "test_plugins_actsvg_detector.svg");
}