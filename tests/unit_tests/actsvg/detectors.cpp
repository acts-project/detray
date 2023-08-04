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
#include "detray/plugins/actsvg_visualization/svg.hpp"
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

int main(int, char**) {

    // The axis
    actsvg::style::stroke stroke_black = actsvg::style::stroke();
    // x-y axis
    auto xy_axis =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250}, stroke_black, "x", "y");  
    // z-r axis
    auto zr_axis =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250}, stroke_black, "z", "r");  

    // Creating the detector and geomentry context.
    using toy_detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    toy_detector_t::geometry_context context{};

    // Creating the converter for the detector.
    detray::actsvg_visualization::svg_converter converter{det, names};

    // Get the svg of the toy detetector in x-y view.
    const auto svg_xy = converter.xy_detector("detector_xy", context);
    // Write the svg of toy detector.
    detray::actsvg_visualization::write_svg("test_actsvg_detector_xy.svg", {xy_axis, svg_xy});

    // Get the svg of the toy detetector in z-r view.
    const auto svg_zr = converter.zr_detector("detector_zr", context);
    // Write the svg of toy detector in z-r view
    detray::actsvg_visualization::write_svg("test_actsvg_detector_zr.svg", {zr_axis, svg_zr});

}