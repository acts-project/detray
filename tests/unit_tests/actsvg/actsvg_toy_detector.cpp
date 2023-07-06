
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
    // The 13th surface in the detector should be a disc
    const auto disc_descr = det.surface_lookup()[13];
    const auto disc_surface = detray::surface{det, disc_descr};
    auto proto_surface = disc_surface.visit_mask<to_proto_surface<>>();

    // Style proto surface.
    proto_surface._fill = style::fill({{0, 100, 0}, 0.5});

    // Draw proto surface.
    views::x_y x_y_view;
    const auto svg = display::surface("surface", proto_surface, x_y_view);

    // Draw x-y-axis.
    style::stroke stroke_black = style::stroke();
    auto x_y_a = draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");

    // Create SVG file.
    svg::file file;
    file.add_object(svg);
    file.add_object(x_y_a);

    // Write SVG File.
    std::ofstream stream;
    stream.open("detray_actsvg.svg");
    stream << file;
    stream.close();
}