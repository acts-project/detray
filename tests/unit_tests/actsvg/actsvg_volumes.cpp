// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/actsvg/svg_conversion.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/core/defs.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/meta.hpp"

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

    // Draw x-y-axis.
    style::stroke stroke_black = style::stroke();
    auto axis =
        draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");  


    for (size_t i = 0; i < det.volumes().size(); i++){
        const auto name = std::string("test_plugins_actsvg_volume") + std::to_string(i);
        // Draw volume
        const auto svg = detray::actsvg_visualization::svg(name, det, det.volume_by_index(i), context, view);
        detray::actsvg_visualization::write_svg(std::array{axis, svg}, name + ".svg");
    }
}