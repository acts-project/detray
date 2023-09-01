/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/utils/groups.hpp"
#include "detray/plugins/svgtools/writer.hpp"
#include "tests/common/tools/particle_gun.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/web/web_builder.hpp"

// System include(s)
#include <array>
#include <filesystem>
#include <string>

int main(int, char**) {

    // In this test we will create a web page to show the detector geometry and
    // more. We will start by creating the svgs we want to include in the web
    // page.

    // Axes.
    const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                             actsvg::style::stroke());

    // Creating the views.
    const actsvg::views::x_y view;

    // Creating the detector and geomentry context.
    using detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    detector_t::geometry_context context{};

    using transform3_t = typename detector_t::transform3;
    using vector3 = typename detector_t::vector3;

    // Creating the svg generator for the detector.
    const detray::svgtools::illustrator il{det, context, true};

    // The vector of svgs that we want to include on the webpage.
    std::vector<actsvg::svg::object> svgs;

    // Indexes of the volumes in the detector to be visualized.
    std::array indices{0UL,  1UL,  2UL,  3UL,  4UL,  5UL,  6UL,
                       7UL,  8UL,  9UL,  10UL, 11UL, 12UL, 13UL,
                       14UL, 15UL, 16UL, 17UL, 18UL, 19UL};

    // Draw the volumes and include them in the svg vector.
    for (std::size_t i : indices) {
        std::string name = "Volume " + std::to_string(i);
        const auto svg = il.draw_volume(name, i, view);
        svgs.push_back(svg);
    }

    // Draw the grids and include them in the svg vector.
    for (const auto i : indices) {
        if (i % 2 == 0) {
            continue;
        }
        std::string name = "Grid " + std::to_string(i);
        const auto svg = il.draw_grid(name, i, view);
        svgs.push_back(svg);
    }

    // Draw some example trajectories and include them in the svg vector (along
    // with their intersections).
    for (const auto qop : std::vector{-4, -8, -16}) {
        std::string name = "Helix (qop: " + std::to_string(qop) + ")";

        const typename detector_t::point3 ori{0.f, 0.f, 80.f};
        const typename detector_t::point3 dir{0, 1, 1};

        // Create the helix trajectory.
        // Constant magnetic field
        vector3 B{0.f * detray::unit<detray::scalar>::T,
                  0.f * detray::unit<detray::scalar>::T,
                  1.f * detray::unit<detray::scalar>::T};

        const detray::detail::helix<transform3_t> helix(
            ori, 0.f, dir, static_cast<float>(qop), &B);
        const auto helix_ir = detray::particle_gun::shoot_particle(det, helix);

        // Draw the helix trajectory.
        const auto svg_helix =
            il.draw_trajectory(name + "_trajectory", helix, view);

        // Draw the intersection record.
        const auto svg_helix_ir =
            il.draw_intersections(name + "_record", helix_ir, view);

        // We one the trajectory and intersection record to be considered as one
        // svg. Thus we group them together before adding the group to the svg
        // vector.
        const auto svg_group = detray::svgtools::utils::group(
            name, std::vector{svg_helix, svg_helix_ir});
        svgs.push_back(svg_group);
    }

    // The output directory for the web page.
    const auto current_directory = std::filesystem::current_path();

    // Create the web page builder.
    actsvg::web::web_builder builder{};

    // To visualize the svg objects in a specific order we need to pass a
    // comparator before we build the page. For instance we might want to
    // display helices on top on the volumes (rather than the other way around).
    // In this example we will alphanumeric_compare which renderes the svgs such
    // that the id of the svg object with the greatest id alphanumerically will
    // be displayed on top. Build the web page.
    auto alphanum_cmp = actsvg::web::compare::alphanumeric;
    builder.build(current_directory / "test_svgtools_website", svgs,
                  alphanum_cmp);

    // Once the direcroy has been created, run the server using "python3 -m
    // http.server" in the directory. Subsequently connect to localhost using
    // the respective port. On the web page, click and drag to move the view
    // box. Use the scroll wheel to zoom.
}
