
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/surface.hpp"
#include "core/detector.hpp"
#include "masks/cylinder3.hpp"
#include "masks/ring2.hpp"
#include "tests/common/test_surfaces.hpp"

#include <string>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

// Three-dimensional definitions
using transform3 = __plugin::transform3;
using vector3 = __plugin::transform3::vector3;
using point3 = __plugin::transform3::point3;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();
constexpr scalar isclose = 1e-5;

using cdetector = detector<transform3>;

cdetector createDetector() {
    cdetector d("three_layer_detector");

    // An inner volume: call it pb
    scalar bp_radius = 29.;
    scalar bp_length = 1000.;
    scalar px_barrel = 600.;
    scalar px_endcap = 0.5 * (bp_length - px_barrel);

    cdetector::volume &bp = d.new_volume("bp", {0., bp_radius, -0.5 * bp_length, 0.5 * bp_length});
    cdetector::portal_cylinder_mask bp_c_ecn = {{bp_radius, -0.5 * bp_length, -0.5 * px_barrel}, {0, 1, -1}};
    cdetector::portal_cylinder_mask bp_c_b = {{bp_radius, -0.5 * px_barrel, 0.5 * px_barrel}, {0, 2, -1}};
    cdetector::portal_cylinder_mask bp_c_ecp = {{bp_radius, 0.5 * px_barrel, 0.5 * bp_length}, {0, 3, -1}};
    cdetector::portal_disc_mask bp_n_disc = {{0., bp_radius}, {-1, 0, -1}};
    cdetector::portal_disc_mask bp_p_disc = {{0., bp_radius}, {0, -1, -1}};
    dvector<cdetector::portal_cylinder_mask> bp_c_portals = {bp_c_ecn, bp_c_b, bp_c_ecp};
    d.add_portal_surface<cdetector::portal_cylinder_mask>(std::move(transform3()), bp_c_portals, bp);
    d.add_portal_surface<cdetector::portal_disc_mask>(std::move(transform3(vector3{0., 0., -0.5 * bp_length})), {bp_n_disc}, bp);
    d.add_portal_surface<cdetector::portal_disc_mask>(std::move(transform3(vector3{0., 0., 0.5 * bp_length})), {bp_p_disc}, bp);
    // Insert an actual beam pipe
    darray<scalar, 3> bpm_values = {25., -0.5 * bp_length + 1., 0.5 * bp_length - 1.};
    d.add_surfaces<cdetector::cylinder_mask>({transform3()}, bpm_values, bp);

    // A wrapping volume, let's call it px
    scalar px_inner_radius = bp_radius;
    scalar px_outer_radius = 55.;

    cdetector::volume &px_ecn = d.new_volume("px_ecn", {px_inner_radius, px_outer_radius, -bp_length, -px_barrel});
    cdetector::portal_cylinder_mask px_ecn_inner = {{px_inner_radius, -bp_length, -px_barrel}, {0, 1, -1}};
    cdetector::portal_cylinder_mask px_ecn_outer = {{px_outer_radius, -bp_length, -px_barrel}, {1, -1, -1}};
    cdetector::portal_disc_mask px_ecn_ecn = {{px_inner_radius, px_outer_radius}, {-1, 1, -1}};
    cdetector::portal_disc_mask px_ecn_ecp = {{px_inner_radius, px_outer_radius}, {1, 2, -1}};
    d.add_portal_surface<cdetector::portal_cylinder_mask>(std::move(transform3()), {px_ecn_inner, px_ecn_outer}, px_ecn);
    d.add_portal_surface<cdetector::portal_disc_mask>(std::move(transform3(vector3{0., 0., -0.5 * bp_length})), {px_ecn_ecn}, px_ecn);
    d.add_portal_surface<cdetector::portal_disc_mask>(std::move(transform3(vector3{0., 0., -0.5 * px_barrel})), {px_ecn_ecp}, px_ecn);
    auto endcap_n = endcap_description(29, 50, -px_barrel - 25., 1., 6, 0.2);
    d.add_surfaces<cdetector::trapezoid_mask>(std::get<dvector<transform3>>(endcap_n),
                                              std::get<cdetector::trapezoid_mask::mask_values>(endcap_n),
                                              px_ecn);

    cdetector::volume &px_b = d.new_volume("px_b", {px_inner_radius, px_outer_radius, -px_barrel, px_barrel});
    cdetector::portal_cylinder_mask px_b_inner = {{px_inner_radius, -px_barrel, px_barrel}, {0, 2, -1}};
    cdetector::portal_cylinder_mask px_b_outer = {{px_outer_radius, -px_barrel, px_barrel}, {2, -1, -1}};
    cdetector::portal_disc_mask px_b_ecn = {{px_inner_radius, px_outer_radius}, {1, 2, -1}};
    cdetector::portal_disc_mask px_b_ecp = {{px_inner_radius, px_outer_radius}, {2, 3, -1}};
    d.add_portal_surface<cdetector::portal_cylinder_mask>(std::move(transform3()), {px_b_inner, px_b_outer}, px_b);
    d.add_portal_surface<cdetector::portal_disc_mask>(std::move(transform3(vector3{0., 0., -0.5 * px_barrel})), {px_b_ecn}, px_b);
    d.add_portal_surface<cdetector::portal_disc_mask>(std::move(transform3(vector3{0., 0., +0.5 * px_barrel})), {px_b_ecp}, px_b);
    auto barrel = barrel_description(33, 0.25, 12, 0.12, 0.25, 600, 2., 7);
    d.add_surfaces<cdetector::rectangle_mask>(std::get<dvector<transform3>>(barrel),
                                              std::get<cdetector::rectangle_mask::mask_values>(barrel),
                                              px_b);

    cdetector::volume &px_ecp = d.new_volume("px_ecp", {px_inner_radius, px_outer_radius, px_barrel, bp_length});
    cdetector::portal_cylinder_mask px_ecp_inner = {{px_inner_radius, px_barrel, bp_length}, {0, 3, -1}};
    cdetector::portal_cylinder_mask px_ecp_outer = {{px_outer_radius, px_barrel, bp_length}, {3, -1, -1}};
    cdetector::portal_disc_mask px_ecp_ecn = {{px_inner_radius, px_outer_radius}, {2, 3, -1}};
    cdetector::portal_disc_mask px_ecp_ecp = {{px_inner_radius, px_outer_radius}, {3, -1, -1}};
    d.add_portal_surface<cdetector::portal_cylinder_mask>(std::move(transform3()), {px_ecp_inner, px_ecp_outer}, px_ecp);
    d.add_portal_surface<cdetector::portal_disc_mask>(std::move(transform3(vector3{0., 0., 0.5 * px_barrel})), {px_ecp_ecn}, px_ecp);
    d.add_portal_surface<cdetector::portal_disc_mask>(std::move(transform3(vector3{0., 0., 0.5 * bp_length})), {px_ecp_ecp}, px_ecp);
    auto endcap_p = endcap_description(29, 50, px_barrel + 25., 1., 6, 0.2);
    d.add_surfaces<cdetector::trapezoid_mask>(std::get<dvector<transform3>>(endcap_p),
                                              std::get<cdetector::trapezoid_mask::mask_values>(endcap_p),
                                              px_ecp);

    return d;
}