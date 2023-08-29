/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <array>
#include <string>

#include "detray/intersection/detail/trajectories.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/particle_gun.hpp"

int main(int, char**) {

    // Axes.
    const auto axes =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                               actsvg::style::stroke(), "axis1", "axis2");

    // Creating the view.
    const actsvg::views::z_r view;

    // Creating the detector and geomentry context.
    using detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    detector_t::geometry_context context{};

    using transform3_t = typename detector_t::transform3;

    // Svg for the detector.
    const detray::svgtools::illustrator il{det, context};
    const auto svg_det = il.draw_detector("detector", view);

    // Creating the rays.
    unsigned int theta_steps{10u};
    unsigned int phi_steps{10u};
    const typename detector_t::point3 ori{0.f, 0.f, 100.f};

    size_t index = 0;
    // Iterate through uniformly distributed momentum directions with ray
    for (const auto test_ray :
         detray::uniform_track_generator<detray::detail::ray<transform3_t>>(
             theta_steps, phi_steps, ori)) {

        // Record all intersections and objects along the ray
        const auto intersection_record =
            detray::particle_gun::shoot_particle(det, test_ray);

        const std::string name = "test_svgtools_intersection_record" + std::to_string(index);
        const auto svg_ir =
            il.draw_intersections(name, intersection_record, view);

        detray::svgtools::write_svg(name + ".svg", {axes, svg_det, svg_ir});

        index++;
    }
}