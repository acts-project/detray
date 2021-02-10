
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

    cdetector d("single_layer_detector");

    scalar barrel_inner_r = 26.;
    scalar barrel_outer_r = 40.;
    scalar barrel_half_length = 600.;

    auto beampipe_index = d.new_volume("beampipe", {0., barrel_inner_r, -barrel_half_length, barrel_half_length});
    cdetector::portal_links beampipe_outer_links = { (beampipe_index), dindex_invalid, dindex_invalid, dindex_invalid};
    cdetector::portal_cylinder_mask beampipe_outer = {{barrel_inner_r, -barrel_half_length, barrel_half_length}, beampipe_outer_links};
    cdetector::portal_links beampipe_ecn_links = {dindex_invalid, (beampipe_index), dindex_invalid, dindex_invalid};
    cdetector::portal_disc_mask beampipe_ecn = {{0., barrel_inner_r}, beampipe_ecn_links};
    cdetector::portal_links beampipe_ecp_links = {(beampipe_index), dindex_invalid, dindex_invalid, dindex_invalid};
    cdetector::portal_disc_mask beampipe_ecp = {{0., barrel_inner_r}, beampipe_ecp_links};

    darray<dindex, 3> beampipe_portal_indices = { dindex_invalid, dindex_invalid, dindex_invalid };
    beampipe_portal_indices[0] = d.add_portal_surface<cdetector::portal_cylinder_mask>(beampipe_index, std::move(transform3()), {beampipe_outer});
    beampipe_portal_indices[1] = d.add_portal_surface<cdetector::portal_disc_mask>(beampipe_index, std::move(transform3(vector3{0., 0., -barrel_half_length})), {beampipe_ecn});
    beampipe_portal_indices[2] = d.add_portal_surface<cdetector::portal_disc_mask>(beampipe_index, std::move(transform3(vector3{0., 0., +barrel_half_length})), {beampipe_ecp});

    auto barrel_index = d.new_volume("single_layer_barrel", {barrel_inner_r, barrel_outer_r, -barrel_half_length, barrel_half_length});
    auto barrel_components = create_barrel_components(32, 0.25, 12, 0.12, 0.25, 2*588, 2., 7, barrel_inner_r, barrel_outer_r, barrel_half_length);
    auto barrel_module = std::get<0>(barrel_components);
    auto barrel_transforms = std::get<1>(barrel_components);
    auto barrel_finders = std::get<2>(barrel_components);

    auto barrel_finder_range = d.add_surface_finders(barrel_index, std::move(barrel_finders));

    cdetector::portal_links barrel_inner_links = { dindex_invalid, (barrel_index), dindex_invalid, barrel_finder_range[0]};
    d.reuse_portal_surface(barrel_index, beampipe_portal_indices[0], barrel_inner_links);

    cdetector::portal_links barrel_outer_links = { (barrel_index), dindex_invalid, barrel_finder_range[1], dindex_invalid};
    cdetector::portal_cylinder_mask barrel_outer = {{barrel_outer_r, -barrel_half_length, barrel_half_length}, beampipe_outer_links};
    cdetector::portal_links barrel_ecn_links = {dindex_invalid, (barrel_index), dindex_invalid, barrel_finder_range[2]};
    cdetector::portal_disc_mask barrel_ecn = {{barrel_inner_r, barrel_outer_r}, barrel_ecn_links};
    cdetector::portal_links barrel_ecp_links = {(barrel_index), dindex_invalid, barrel_finder_range[3], dindex_invalid};
    cdetector::portal_disc_mask barrel_ecp = {{barrel_inner_r, barrel_outer_r}, barrel_ecp_links};

    darray<dindex, 4> barrel_portal_indices = { beampipe_portal_indices[0], dindex_invalid, dindex_invalid, dindex_invalid };
    barrel_portal_indices[1] = d.add_portal_surface<cdetector::portal_cylinder_mask>(barrel_index, std::move(transform3()), {barrel_outer});
    barrel_portal_indices[2] = d.add_portal_surface<cdetector::portal_disc_mask>(barrel_index, std::move(transform3(vector3{0., 0., -barrel_half_length})), {barrel_ecn});
    barrel_portal_indices[3] = d.add_portal_surface<cdetector::portal_disc_mask>(barrel_index, std::move(transform3(vector3{0., 0., +barrel_half_length})), {barrel_ecp});

    return d;

};

