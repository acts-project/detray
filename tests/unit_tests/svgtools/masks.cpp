/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"

// System include(s)
#include <array>
#include <string>

int main(int, char**) {

    // Axes.
    const auto axes =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                               actsvg::style::stroke(), "axis1", "axis2");

    using point3 = std::array<actsvg::scalar, 3>;
    using point3_container = std::vector<point3>;

    using toy_detector_t = detray::detector<detray::toy_metadata<>>;

    const typename toy_detector_t::transform3 transform;
    const actsvg::views::x_y view{};

    // e_min_r, e_max_r, e_min_phi_rel, e_max_phi_rel, e_average_phi, e_shift_x,
    // e_shift_y, e_size
    detray::mask<detray::annulus2D<>> ann2D{0u,  100.0, 200.0, -1.0,
                                            1.0, 0.0,   0.0,   100.0};
    const auto ann2D_proto =
        detray::svgtools::conversion::surface<point3_container>(transform,
                                                                ann2D);
    const auto ann2D_svg = actsvg::display::surface("", ann2D_proto, view);
    detray::svgtools::write_svg("svgtools_annulus2D.svg", {axes, ann2D_svg});

    // e_r, e_n_half_z, e_p_half_z, e_size
    detray::mask<detray::cylinder2D<>> cyl2D{0u, 100.0, -10.0, 10.0};
    const auto cyl2D_proto =
        detray::svgtools::conversion::surface<point3_container>(transform,
                                                                cyl2D);
    const auto cyl2D_svg = actsvg::display::surface("", cyl2D_proto, view);
    detray::svgtools::write_svg("svgtools_cylinder2D.svg", {axes, cyl2D_svg});

    // e_half_x, e_half_y, e_size
    detray::mask<detray::rectangle2D<>> rec2D{0u, 100.0, 100.0};
    const auto rec2D_proto =
        detray::svgtools::conversion::surface<point3_container>(transform,
                                                                rec2D);
    const auto rec2D_svg = actsvg::display::surface("", rec2D_proto, view);
    detray::svgtools::write_svg("svgtools_rectangle2D.svg", {axes, rec2D_svg});

    // e_inner_r, e_outer_r, e_size
    detray::mask<detray::ring2D<>> rin2D{0u, 50.0, 100.0};
    const auto rin2D_proto =
        detray::svgtools::conversion::surface<point3_container>(transform,
                                                                rin2D);
    const auto rin2D_svg = actsvg::display::surface("", rin2D_proto, view);
    detray::svgtools::write_svg("svgtools_ring2D.svg", {axes, rin2D_svg});

    // e_half_length_0, e_half_length_1, e_half_length_2, e_divisor, e_size
    detray::mask<detray::trapezoid2D<>> tra2D{0u, 100.0, 50.0, 200.0};
    const auto tra2D_proto =
        detray::svgtools::conversion::surface<point3_container>(transform,
                                                                tra2D);
    const auto tra2D_svg = actsvg::display::surface("", tra2D_proto, view);
    detray::svgtools::write_svg("svgtools_trapezoid2D.svg", {axes, tra2D_svg});
}