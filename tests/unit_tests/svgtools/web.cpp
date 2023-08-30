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
#include "detray/intersection/detail/trajectories.hpp"
#include "tests/common/tools/particle_gun.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <array>
#include <string>

int main(int, char**) {

    // Axes.
    const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                             actsvg::style::stroke());

    // Creating the views.
    const actsvg::views::z_r view;

    // Creating the detector and geomentry context.
    using detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    detector_t::geometry_context context{};

    using transform3_t = typename detector_t::transform3;
    using vector3 = typename detector_t::vector3;

    // Creating the svg generator for the detector.
    const detray::svgtools::illustrator il{det, context, true};

    // Indexes of the volumes in the detector to be visualized.
    std::array indices{0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL, 8UL, 9UL, 10UL, 11UL, 12UL, 13UL, 14UL, 15UL, 16UL, 17UL, 18UL, 19UL};
    for (std::size_t i : indices) {
        std::string name = "Volume " + std::to_string(i);
        // Visualization of volume i:
        const auto svg = il.draw_volume("", i, view);
        detray::svgtools::write_svg(name + ".svg", svg);
    }

    for (const auto i : indices){
        if (i % 2 == 0){
            continue;
        }
        std::string name = "Grid " + std::to_string(i);
        const auto svg = il.draw_grid(name, i, view);
        detray::svgtools::write_svg(name + ".svg", svg);
    }

    for (const auto qop : std::vector{-4, -8, -16}){
        std::string name = "Helix (qop: " + std::to_string(qop) + ")";

        const typename detector_t::point3 ori{0.f, 0.f, 80.f};
        const typename detector_t::point3 dir{0, 1, 1};

        // Constant magnetic field
        vector3 B{0.f * detray::unit<detray::scalar>::T,
                0.f * detray::unit<detray::scalar>::T,
                1.f * detray::unit<detray::scalar>::T};

        const detray::detail::helix<transform3_t> helix(ori, 0.f, dir, static_cast<float>(qop), &B);
        const auto helix_ir = detray::particle_gun::shoot_particle(det, helix);

        const auto svg_helix = il.draw_trajectory("trajectory", helix, view);
        const auto svg_helix_ir =
            il.draw_intersections("record", helix_ir, view);
        detray::svgtools::write_svg(name + ".svg",
                                    {svg_helix, svg_helix_ir});

    }
}