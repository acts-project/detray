/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"
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
    const auto axes =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250}, actsvg::style::stroke(), "axis1", "axis2");

    // Creating the view.
    const actsvg::views::z_r view;

    // Creating the detector and geomentry context.
    using detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    detector_t::geometry_context context{};

    // Svg for the detector.
    const detray::svgtools::illustrator il{det, names};
    const auto svg_det = il.draw_detector("detector", context, view);
    
    // Creating the rays.
    using transform3_t = typename detector_t::transform3;
    using vector3 =  typename detector_t::vector3;

    const typename detector_t::point3 ori{0.f, 0.f, 100.f};
    const typename detector_t::point3 dir{0,1,1};

    
    const detray::detail::ray<transform3_t> ray(ori, 0.f, dir, 0.f);
    const auto ray_ir =
            detray::particle_gun::shoot_particle(det, ray);
            
    const auto svg_ray = il.draw_trajectory("trajectory", ray, view);
    const auto svg_ray_ir = il.draw_intersections("record", context, ray_ir, view);
    detray::svgtools::write_svg("ray.svg", {svg_det, svg_ray, svg_ray_ir});
    
    // Constant magnetic field
    vector3 B{1.f * detray::unit<detray::scalar>::T, 1.f * detray::unit<detray::scalar>::T,
              1.f * detray::unit<detray::scalar>::T};

    const detray::detail::helix<transform3_t> helix(ori, 0.f, dir, -500.f, &B);
        const auto helix_ir =
            detray::particle_gun::shoot_particle(det, helix);
    
    const auto svg_helix = il.draw_trajectory("trajectory", helix, view);
    const auto svg_helix_ir = il.draw_intersections("record", context, helix_ir, view);
    detray::svgtools::write_svg("helix.svg", {svg_det, svg_helix, svg_helix_ir});

}