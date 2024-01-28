/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/navigation/detail/trajectories.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/utils/particle_gun.hpp"
#include "detray/tracks/tracks.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <string>

GTEST_TEST(svgtools, intersections) {

    // This test creates the visualization using the illustrator class.
    // However, for full control over the process, it is also possible to use
    // the tools in svgstools::conversion, svgstools::display, and
    // actsvg::display by converting the object to a proto object, optionally
    // styling it, and then displaying it.

    // Axes.
    const auto axes =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                               actsvg::style::stroke(), "axis1", "axis2");

    // Creating the view.
    const actsvg::views::z_r view;

    // Creating the detector and geomentry context.
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr);
    using detector_t = decltype(det);

    using transform3_t = typename detector_t::transform3;

    // Creating the illustrator.
    const detray::svgtools::illustrator il{det, names};

    // Drawing the detector.
    const auto svg_det = il.draw_detector(view);

    // Creating the rays.
    using generator_t =
        detray::uniform_track_generator<detray::detail::ray<transform3_t>>;
    auto trk_gen_cfg = generator_t::configuration{};
    trk_gen_cfg.origin({0.f, 0.f, 100.f}).phi_steps(10u).theta_steps(10u);

    std::size_t index = 0;
    // Iterate through uniformly distributed momentum directions with ray
    for (const auto test_ray : generator_t{trk_gen_cfg}) {

        // Record all intersections and objects along the ray
        const auto intersection_record =
            detray::particle_gun::shoot_particle(det, test_ray);

        const std::string name =
            "test_svgtools_intersection_record" + std::to_string(index);

        // Drawing the intersections.
        const auto svg_ir = il.draw_intersections(name, intersection_record,
                                                  test_ray.dir(), view);

        detray::svgtools::write_svg(name, {axes, svg_det, svg_ir});

        index++;
    }
}
