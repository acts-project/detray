
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "actsvg/core.hpp"
#include "actsvg/core/defs.hpp"

// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detector_writer.hpp"


// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/masks/cylinder2D.hpp"
#include "detray/masks/masks.hpp"
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracks/tracks.hpp"
#include "actsvg/meta.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/proto/surface.hpp"

#include <type_traits>

#include "./surface_svg_converter.hpp"


int main(int argc, char* argv[]) {

    using toy_detector_t = detray::detector<detray::toy_metadata<>>;

    vecmem::host_memory_resource host_mr;
    const toy_detector_t det = detray::create_toy_geometry(host_mr, 4, 3);
    detray::surface<detray::detector<detray::toy_metadata<>>>::context context{};
    
    views::x_y x_y_view;
    // Draw x-y-axis.
    style::stroke stroke_black = style::stroke();
    auto x_y_a = draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");

    // Create SVG file.
    svg::file file;

    file.add_object(x_y_a);

    auto indexes = {13, 20, 100, 150, 200, 250};

    for (int i : indexes)
    {
        const auto description = det.surface_lookup()[i];
        const auto surface = detray::surface{det, description};
        auto p_surface = surface_converter::convert_surface(surface, context);

        // Style proto surface.
        p_surface._fill = style::fill({{0, 100, 0}, 0.5});
        // Draw proto surface.
        auto name = "detector_surface" + std::to_string(i);
        const auto svg = display::surface(name, p_surface, x_y_view);

        file.add_object(svg);
    }

    // Write SVG File.
    std::ofstream stream;
    stream.open("test_plugins_actsvg_detector.svg");
    stream << file;
    stream.close();
}