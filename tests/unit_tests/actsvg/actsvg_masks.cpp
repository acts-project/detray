
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
#include "actsvg/styles/defaults.hpp"

#include <type_traits>

#include "./surface_svg_converter.hpp"

using point3 = std::array<scalar, 3>;
using point3_container = std::vector<point3>;
using proto_surface = proto::surface<point3_container>;

void write_svg(proto_surface p_surface, std::string file_name){
    // Style proto surface.
    p_surface._fill = style::fill({{0, 100, 0}, 0.5});

    // Draw proto surface.
    views::x_y x_y_view;
    const auto svg = display::surface("surface", p_surface, x_y_view);

    // Draw x-y-axis.
    style::stroke stroke_black = style::stroke();
    auto x_y_a = draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");

    // Create SVG file.
    svg::file file;
    file.add_object(svg);
    file.add_object(x_y_a);

    // Write SVG File.
    std::ofstream stream;
    stream.open(file_name);
    stream << file;
    stream.close();
}

int main(int argc, char* argv[]) {

    //e_min_r, e_max_r, e_min_phi_rel, e_max_phi_rel, e_average_phi, e_shift_x, e_shift_y, e_size
    detray::mask<detray::annulus2D<>> ann2D{0u, 100.0, 200.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    write_svg(surface_converter::convert_mask(ann2D), "test_plugins_actsvg_annulus2D.svg");

    //e_r, e_n_half_z, e_p_half_z, e_size
    detray::mask<detray::cylinder2D<>> cyl2D{0u, 100.0, -10.0, 10.0};
    write_svg(surface_converter::convert_mask(cyl2D), "test_plugins_actsvg_cylinder2D.svg");

    //e_half_x, e_half_y, e_size
    detray::mask<detray::rectangle2D<>> rec2D{0u, 100.0, 100.0};
    auto p = surface_converter::convert_mask(rec2D);
    auto tf = actsvg::style::transform();
    tf._rot = {90, 100., 0.8};
    tf._tr = { 10, 20};
    p._transform = tf;
    write_svg(p, "test_plugins_actsvg_rectangle2D.svg");

    //e_inner_r, e_outer_r, e_size
    detray::mask<detray::ring2D<>> rin2D{0u, 50.0, 100.0};
    write_svg(surface_converter::convert_mask(rin2D), "test_plugins_actsvg_ring2D.svg");

    //e_half_length_0, e_half_length_1, e_half_length_2, e_divisor, e_size
    detray::mask<detray::trapezoid2D<>> tra2D{0u, 100.0, 50.0, 200.0 };
    write_svg(surface_converter::convert_mask(tra2D), "test_plugins_actsvg_trapezoid2D.svg");

    //e_min_x, e_min_y, e_min_z, e_max_x, e_max_y, e_max_z, e_size
    //detray::mask<detray::cuboid3D<>> cub3D{0u, 0.0, 0.0, 100.0, 100.0, 100.0};
    //write_svg(surface_converter::convert_mask(cub3D), "test_plugins_actsvg_cuboid3D.svg");

    //e_min_r, e_min_phi, e_min_z, e_max_r, e_max_phi, e_max_z, e_size
    //detray::mask<detray::cylinder3D> cyl3D{0u, 100.0, 0.0, 0.0, 300.0, 4.0, 10.0};
    //write_svg(surface_converter::convert_mask(cyl3D), "test_plugins_actsvg_cylinder3D.svg");
}
