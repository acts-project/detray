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
#include "detray/plugins/actsvg/surface_conversion.hpp"
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

    views::x_y x_y_view;

    // Draw x-y-axis.
    style::stroke stroke_black = style::stroke();
    auto x_y_a =
        draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");

    // Create SVG file.
    svg::file file;

    file.add_object(x_y_a);

    auto indices = {13, 20, 100, 150, 200, 250};

    for (const auto& pair :
         detray::views::pick(det.surface_lookup(), indices)) {
        const auto index = pair.first;
        const auto description = pair.second;

        const auto surface = detray::surface{det, description};
        auto p_surface =
            detray::actsvg_visualization::convert_surface(surface, context);

        // Style proto surface.
        p_surface._fill = style::fill({{0, 100, 0}, 0.5});

        // Draw proto surface.
        auto name = "detector_surface" + std::to_string(index);
        const auto svg = display::surface(name, p_surface, x_y_view);

        file.add_object(svg);
    }

    // Write SVG File.
    detray::io::detail::file_handle stream{"test_plugins_actsvg_detector",
                                           ".svg",
                                           std::ios::out | std::ios::trunc};
    *stream << file;
}