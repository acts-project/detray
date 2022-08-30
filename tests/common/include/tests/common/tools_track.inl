/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/transform_store.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/masks/rectangle2.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/tracks/tracks.hpp"

// Google Test include(s)
#include <gtest/gtest.h>

using namespace detray;
using vector2 = __plugin::vector2<scalar>;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using transform3 = __plugin::transform3<detray::scalar>;
using matrix_operator = standard_matrix_operator<scalar>;
using track_helper = detail::track_helper<matrix_operator>;

enum mask_ids : unsigned int {
    e_rectangle2 = 0,
};

enum material_ids : unsigned int {
    e_slab = 0,
};

using mask_defs = tuple_vector_registry<mask_ids, rectangle2<>>;
using material_defs =
    tuple_vector_registry<material_ids, material_slab<scalar>>;

TEST(tools, bound_track_parameters) {

    // surface container
    std::vector<surface<mask_defs, material_defs>> surfaces;
    surfaces.emplace_back(0, mask_defs::link_type{e_rectangle2, 0},
                          material_defs::link_type{e_slab, 0}, 0, 0, false);
    surfaces.emplace_back(1, mask_defs::link_type{e_rectangle2, 0},
                          material_defs::link_type{e_slab, 0}, 0, 0, false);

    /// Declare track parameters

    // first track
    dindex sf_idx1 = 0;
    typename bound_track_parameters<transform3>::vector_type param1;
    getter::element(param1, e_bound_loc0, 0) = 1.;
    getter::element(param1, e_bound_loc1, 0) = 2.;
    getter::element(param1, e_bound_phi, 0) = 0.1;
    getter::element(param1, e_bound_theta, 0) = 0.2;
    getter::element(param1, e_bound_qoverp, 0) = -0.01;
    getter::element(param1, e_bound_time, 0) = 0.1;

    typename bound_track_parameters<transform3>::covariance_type cov1 =
        matrix_operator().template zero<e_bound_size, e_bound_size>();

    bound_track_parameters<transform3> trck1(sf_idx1, param1, cov1);

    // second track
    dindex sf_idx2 = 1;
    typename bound_track_parameters<transform3>::vector_type param2;
    getter::element(param2, e_bound_loc0, 0) = 4.;
    getter::element(param2, e_bound_loc1, 0) = 20.;
    getter::element(param2, e_bound_phi, 0) = 0.8;
    getter::element(param2, e_bound_theta, 0) = 1.4;
    getter::element(param2, e_bound_qoverp, 0) = 1.;
    getter::element(param2, e_bound_time, 0) = 0.;

    typename bound_track_parameters<transform3>::covariance_type cov2 =
        matrix_operator().template zero<e_bound_size, e_bound_size>();

    bound_track_parameters<transform3> trck2(sf_idx2, param2, cov2);

    /// Check the elements

    // first track
    EXPECT_FLOAT_EQ(trck1.local()[0], getter::element(param1, e_bound_loc0, 0));
    EXPECT_FLOAT_EQ(trck1.local()[1], getter::element(param1, e_bound_loc1, 0));
    EXPECT_FLOAT_EQ(trck1.phi(), getter::element(param1, e_bound_phi, 0));
    EXPECT_FLOAT_EQ(trck1.theta(), getter::element(param1, e_bound_theta, 0));
    EXPECT_FLOAT_EQ(trck1.qop(), getter::element(param1, e_bound_qoverp, 0));
    EXPECT_FLOAT_EQ(trck1.charge(), -1);
    EXPECT_FLOAT_EQ(trck1.time(), getter::element(param1, e_bound_time, 0));

    // second track
    EXPECT_FLOAT_EQ(trck2.local()[0], getter::element(param2, e_bound_loc0, 0));
    EXPECT_FLOAT_EQ(trck2.local()[1], getter::element(param2, e_bound_loc1, 0));
    EXPECT_FLOAT_EQ(trck2.phi(), getter::element(param2, e_bound_phi, 0));
    EXPECT_FLOAT_EQ(trck2.theta(), getter::element(param2, e_bound_theta, 0));
    EXPECT_FLOAT_EQ(trck2.qop(), getter::element(param2, e_bound_qoverp, 0));
    EXPECT_FLOAT_EQ(trck2.charge(), 1.);
    EXPECT_FLOAT_EQ(trck2.time(), getter::element(param2, e_bound_time, 0));
}

TEST(tools, free_track_parameters) {

    point3 pos = {4., 10., 2.};
    scalar time = 0.1;
    vector3 mom = {10., 20., 30.};
    scalar charge = -1.;

    typename free_track_parameters<transform3>::vector_type param;
    getter::element(param, e_free_pos0, 0) = pos[0];
    getter::element(param, e_free_pos1, 0) = pos[1];
    getter::element(param, e_free_pos2, 0) = pos[2];
    getter::element(param, e_free_time, 0) = time;
    getter::element(param, e_free_dir0, 0) = mom[0] / getter::norm(mom);
    getter::element(param, e_free_dir1, 0) = mom[1] / getter::norm(mom);
    getter::element(param, e_free_dir2, 0) = mom[2] / getter::norm(mom);
    getter::element(param, e_free_qoverp, 0) = charge / getter::norm(mom);

    typename free_track_parameters<transform3>::covariance_type cov;

    // first constructor
    free_track_parameters<transform3> trck1(param, cov);
    EXPECT_FLOAT_EQ(trck1.pos()[0], getter::element(param, e_free_pos0, 0));
    EXPECT_FLOAT_EQ(trck1.pos()[1], getter::element(param, e_free_pos1, 0));
    EXPECT_FLOAT_EQ(trck1.pos()[2], getter::element(param, e_free_pos2, 0));
    EXPECT_FLOAT_EQ(trck1.dir()[0], getter::element(param, e_free_dir0, 0));
    EXPECT_FLOAT_EQ(trck1.dir()[1], getter::element(param, e_free_dir1, 0));
    EXPECT_FLOAT_EQ(trck1.dir()[2], getter::element(param, e_free_dir2, 0));
    EXPECT_FLOAT_EQ(getter::norm(trck1.mom()), getter::norm(mom));
    EXPECT_FLOAT_EQ(trck1.time(), getter::element(param, e_free_time, 0));
    EXPECT_FLOAT_EQ(trck1.qop(), getter::element(param, e_free_qoverp, 0));
    EXPECT_FLOAT_EQ(trck1.pT(),
                    std::sqrt(std::pow(mom[0], 2) + std::pow(mom[1], 2)));

    // second constructor
    free_track_parameters<transform3> trck2(pos, time, mom, charge);
    EXPECT_FLOAT_EQ(trck2.pos()[0], pos[0]);
    EXPECT_FLOAT_EQ(trck2.pos()[1], pos[1]);
    EXPECT_FLOAT_EQ(trck2.pos()[2], pos[2]);
    EXPECT_FLOAT_EQ(trck2.dir()[0], mom[0] / getter::norm(mom));
    EXPECT_FLOAT_EQ(trck2.dir()[1], mom[1] / getter::norm(mom));
    EXPECT_FLOAT_EQ(trck2.dir()[2], mom[2] / getter::norm(mom));
    EXPECT_FLOAT_EQ(getter::norm(trck2.mom()), getter::norm(mom));
    EXPECT_FLOAT_EQ(trck2.time(), time);
    EXPECT_FLOAT_EQ(trck2.qop(), charge / getter::norm(mom));
    EXPECT_FLOAT_EQ(trck2.pT(),
                    std::sqrt(std::pow(mom[0], 2) + std::pow(mom[1], 2)));
}
