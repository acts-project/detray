// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/surface.hpp"
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
#include "actsvg/styles/defaults.hpp"

// System include(s)
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

using namespace actsvg;

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;
using proto_surface = proto::surface<point3_container>;

int main(int, char**) {

    views::x_y view;

    style::stroke stroke_black = style::stroke();
    auto axis =
        draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");

    // e_min_r, e_max_r, e_min_phi_rel, e_max_phi_rel, e_average_phi, e_shift_x,
    // e_shift_y, e_size
    detray::mask<detray::annulus2D<>> ann2D{0u,  100.0, 200.0, 0.0,
                                            0.0, 0.0,   0.0,   0.0};
    auto name = "test_plugins_actsvg_annulus2D";
    auto svg = detray::actsvg_visualization::svg(name, ann2D, view);
    detray::actsvg_visualization::write_svg(std::array{svg, axis}, name + std::string(".svg"));

    // e_r, e_n_half_z, e_p_half_z, e_size
    detray::mask<detray::cylinder2D<>> cyl2D{0u, 100.0, -10.0, 10.0};
    name = "test_plugins_actsvg_cylinder2D";
    svg = detray::actsvg_visualization::svg(name, cyl2D, view);
    detray::actsvg_visualization::write_svg(std::array{svg, axis}, name + std::string(".svg"));

    // e_half_x, e_half_y, e_size
    detray::mask<detray::rectangle2D<>> rec2D{0u, 100.0, 100.0};
    name = "test_plugins_actsvg_rectangle2D";
    svg = detray::actsvg_visualization::svg(name, rec2D, view);
    detray::actsvg_visualization::write_svg(std::array{svg, axis}, name + std::string(".svg"));

    // e_inner_r, e_outer_r, e_size
    detray::mask<detray::ring2D<>> rin2D{0u, 50.0, 100.0};
    name = "test_plugins_actsvg_ring2D";
    svg = detray::actsvg_visualization::svg(name, rin2D, view);
    detray::actsvg_visualization::write_svg(std::array{svg, axis}, name + std::string(".svg"));

    // e_half_length_0, e_half_length_1, e_half_length_2, e_divisor, e_size
    detray::mask<detray::trapezoid2D<>> tra2D{0u, 100.0, 50.0, 200.0};
    name = "test_plugins_actsvg_trapezoid2D";
    svg = detray::actsvg_visualization::svg(name, tra2D, view);
    detray::actsvg_visualization::write_svg(std::array{svg, axis}, name + std::string(".svg"));

    // Transform test:
    auto transform = actsvg::style::transform();
    transform._rot = {35., 100., 100.};
    transform._tr = {0, 0};

    auto rect2D_proto = detray::actsvg_visualization::convert_mask(rec2D);
    rect2D_proto._transform = transform;
    name = "test_plugins_actsvg_rectangle2D_transform";
    svg = actsvg::display::surface(name, rect2D_proto, view);
    detray::actsvg_visualization::write_svg(std::array{svg, axis}, name + std::string(".svg"));

    auto rin2D_proto = detray::actsvg_visualization::convert_mask(rin2D);
    rin2D_proto._transform = transform;
    name = "test_plugins_actsvg_ring2D_transform";
    svg = actsvg::display::surface(name, rin2D_proto, view);
    detray::actsvg_visualization::write_svg(std::array{svg, axis}, name + std::string(".svg"));
}
