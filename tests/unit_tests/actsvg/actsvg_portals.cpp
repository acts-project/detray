// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/plugins/actsvg/svg_conversion.hpp"

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

    for (const auto& description : det.portals()) {

        const auto portal = detray::surface{det, description};
        const auto name = "toy_detector_portal";

        const auto svg = detray::actsvg_visualization::svg(name, det, portal, context, view);

        svg::file file;
        file.add_object(svg);
        detray::io::detail::file_handle stream{std::string("test_plugins_actsvg_") + name,
                                           ".svg",
                                           std::ios::out | std::ios::trunc};
        *stream << file;
    }
}