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
#include "detray/masks/masks.hpp"
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

using mask_defs = tuple_vector_registry<mask_ids, mask<rectangle2D<>>>;
using material_defs =
    tuple_vector_registry<material_ids, material_slab<scalar>>;
using mask_link_t = dtyped_index<mask_ids, dindex>;
using material_link_t = dtyped_index<material_ids, dindex>;

TEST(tools, bound_track_parameters) {

    // surface container
    std::vector<surface<mask_link_t, material_link_t>> surfaces;
    surfaces.emplace_back(0, mask_link_t{e_rectangle2, 0},
                          material_link_t{e_slab, 0}, 0, 0,
                          surface_id::e_sensitive);
    surfaces.emplace_back(1, mask_link_t{e_rectangle2, 0},
                          material_link_t{e_slab, 0}, 0, 0,
                          surface_id::e_sensitive);

    /// Declare track parameters

    // first track
    dindex sf_idx1 = 0;
    typename bound_track_parameters<transform3>::vector_type bound_vec1 =
        matrix_operator().template zero<e_bound_size, 1>();
    getter::element(bound_vec1, e_bound_loc0, 0) = 1.;
    getter::element(bound_vec1, e_bound_loc1, 0) = 2.;
    getter::element(bound_vec1, e_bound_phi, 0) = 0.1;
    getter::element(bound_vec1, e_bound_theta, 0) = 0.2;
    getter::element(bound_vec1, e_bound_qoverp, 0) = -0.01;
    getter::element(bound_vec1, e_bound_time, 0) = 0.1;

    typename bound_track_parameters<transform3>::covariance_type bound_cov1 =
        matrix_operator().template zero<e_bound_size, e_bound_size>();

    bound_track_parameters<transform3> bound_param1(sf_idx1, bound_vec1,
                                                    bound_cov1);

    // second track
    dindex sf_idx2 = 1;
    typename bound_track_parameters<transform3>::vector_type bound_vec2 =
        matrix_operator().template zero<e_bound_size, 1>();
    getter::element(bound_vec2, e_bound_loc0, 0) = 4.;
    getter::element(bound_vec2, e_bound_loc1, 0) = 20.;
    getter::element(bound_vec2, e_bound_phi, 0) = 0.8;
    getter::element(bound_vec2, e_bound_theta, 0) = 1.4;
    getter::element(bound_vec2, e_bound_qoverp, 0) = 1.;
    getter::element(bound_vec2, e_bound_time, 0) = 0.;

    typename bound_track_parameters<transform3>::covariance_type bound_cov2 =
        matrix_operator().template zero<e_bound_size, e_bound_size>();

    bound_track_parameters<transform3> bound_param2(sf_idx2, bound_vec2,
                                                    bound_cov2);
    bound_track_parameters<transform3> bound_param3(sf_idx2, bound_vec2,
                                                    bound_cov2);

    /// Check the elements

    // first track
    EXPECT_FLOAT_EQ(bound_param1.local()[0],
                    getter::element(bound_vec1, e_bound_loc0, 0));
    EXPECT_FLOAT_EQ(bound_param1.local()[1],
                    getter::element(bound_vec1, e_bound_loc1, 0));
    EXPECT_FLOAT_EQ(bound_param1.phi(),
                    getter::element(bound_vec1, e_bound_phi, 0));
    EXPECT_FLOAT_EQ(bound_param1.theta(),
                    getter::element(bound_vec1, e_bound_theta, 0));
    EXPECT_FLOAT_EQ(bound_param1.qop(),
                    getter::element(bound_vec1, e_bound_qoverp, 0));
    EXPECT_FLOAT_EQ(bound_param1.charge(), -1);
    EXPECT_FLOAT_EQ(bound_param1.time(),
                    getter::element(bound_vec1, e_bound_time, 0));

    // second track
    EXPECT_FLOAT_EQ(bound_param2.local()[0],
                    getter::element(bound_vec2, e_bound_loc0, 0));
    EXPECT_FLOAT_EQ(bound_param2.local()[1],
                    getter::element(bound_vec2, e_bound_loc1, 0));
    EXPECT_FLOAT_EQ(bound_param2.phi(),
                    getter::element(bound_vec2, e_bound_phi, 0));
    EXPECT_FLOAT_EQ(bound_param2.theta(),
                    getter::element(bound_vec2, e_bound_theta, 0));
    EXPECT_FLOAT_EQ(bound_param2.qop(),
                    getter::element(bound_vec2, e_bound_qoverp, 0));
    EXPECT_FLOAT_EQ(bound_param2.charge(), 1.);
    EXPECT_FLOAT_EQ(bound_param2.time(),
                    getter::element(bound_vec2, e_bound_time, 0));

    EXPECT_TRUE(!(bound_param2 == bound_param1));
    EXPECT_TRUE(bound_param2 == bound_param3);
}

TEST(tools, free_track_parameters) {

    point3 pos = {4., 10., 2.};
    scalar time = 0.1;
    vector3 mom = {10., 20., 30.};
    scalar charge = -1.;

    typename free_track_parameters<transform3>::vector_type free_vec =
        matrix_operator().template zero<e_free_size, 1>();
    getter::element(free_vec, e_free_pos0, 0) = pos[0];
    getter::element(free_vec, e_free_pos1, 0) = pos[1];
    getter::element(free_vec, e_free_pos2, 0) = pos[2];
    getter::element(free_vec, e_free_time, 0) = time;
    getter::element(free_vec, e_free_dir0, 0) = mom[0] / getter::norm(mom);
    getter::element(free_vec, e_free_dir1, 0) = mom[1] / getter::norm(mom);
    getter::element(free_vec, e_free_dir2, 0) = mom[2] / getter::norm(mom);
    getter::element(free_vec, e_free_qoverp, 0) = charge / getter::norm(mom);

    typename free_track_parameters<transform3>::covariance_type free_cov =
        matrix_operator().template zero<e_free_size, e_free_size>();

    // first constructor
    free_track_parameters<transform3> free_param1(free_vec, free_cov);
    EXPECT_FLOAT_EQ(free_param1.pos()[0],
                    getter::element(free_vec, e_free_pos0, 0));
    EXPECT_FLOAT_EQ(free_param1.pos()[1],
                    getter::element(free_vec, e_free_pos1, 0));
    EXPECT_FLOAT_EQ(free_param1.pos()[2],
                    getter::element(free_vec, e_free_pos2, 0));
    EXPECT_FLOAT_EQ(free_param1.dir()[0],
                    getter::element(free_vec, e_free_dir0, 0));
    EXPECT_FLOAT_EQ(free_param1.dir()[1],
                    getter::element(free_vec, e_free_dir1, 0));
    EXPECT_FLOAT_EQ(free_param1.dir()[2],
                    getter::element(free_vec, e_free_dir2, 0));
    EXPECT_FLOAT_EQ(getter::norm(free_param1.mom()), getter::norm(mom));
    EXPECT_FLOAT_EQ(free_param1.time(),
                    getter::element(free_vec, e_free_time, 0));
    EXPECT_FLOAT_EQ(free_param1.qop(),
                    getter::element(free_vec, e_free_qoverp, 0));
    EXPECT_FLOAT_EQ(free_param1.pT(),
                    std::sqrt(std::pow(mom[0], 2) + std::pow(mom[1], 2)));

    // second constructor
    free_track_parameters<transform3> free_param2(pos, time, mom, charge);
    EXPECT_FLOAT_EQ(free_param2.pos()[0], pos[0]);
    EXPECT_FLOAT_EQ(free_param2.pos()[1], pos[1]);
    EXPECT_FLOAT_EQ(free_param2.pos()[2], pos[2]);
    EXPECT_FLOAT_EQ(free_param2.dir()[0], mom[0] / getter::norm(mom));
    EXPECT_FLOAT_EQ(free_param2.dir()[1], mom[1] / getter::norm(mom));
    EXPECT_FLOAT_EQ(free_param2.dir()[2], mom[2] / getter::norm(mom));
    EXPECT_FLOAT_EQ(getter::norm(free_param2.mom()), getter::norm(mom));
    EXPECT_FLOAT_EQ(free_param2.time(), time);
    EXPECT_FLOAT_EQ(free_param2.qop(), charge / getter::norm(mom));
    EXPECT_FLOAT_EQ(free_param2.pT(),
                    std::sqrt(std::pow(mom[0], 2) + std::pow(mom[1], 2)));

    EXPECT_TRUE(free_param2 == free_param1);
}
