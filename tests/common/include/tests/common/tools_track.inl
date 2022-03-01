/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/core/transform_store.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/masks/rectangle2.hpp"
#include "detray/tools/track.hpp"

using namespace detray;
using vector2 = __plugin::vector2<detray::scalar>;

enum mask_ids : unsigned int {
    e_rectangle2 = 0,
};

using mask_defs = mask_registry<mask_ids, rectangle2<>>;

TEST(tools, bound_track_parameters) {

    // surface container
    std::vector<surface<mask_defs>> surfaces;
    surfaces.emplace_back(0, mask_defs::link_type{e_rectangle2, 0}, 0, 0,
                          false);
    surfaces.emplace_back(1, mask_defs::link_type{e_rectangle2, 0}, 0, 0,
                          false);

    // transform container
    static_transform_store trfs;
    typename decltype(trfs)::context ctx{};

    vector3 z1 = vector::normalize(vector3{1, 0., 0.});
    vector3 x1 = vector::normalize(vector3{0., 0., 1.});
    vector3 y1 = vector::cross(z1, x1);
    point3 t1{2., 3., 4.};
    trfs.emplace_back(ctx, t1, z1, x1);

    vector3 z2 = vector::normalize(vector3{2., 3., 4.});
    vector3 x2 = vector::normalize(vector3{3., -2., 0.});
    vector3 y2 = vector::cross(z2, x2);
    point3 t2{10., 0., 6.};
    trfs.emplace_back(ctx, t2, z2, x2);

    /// Declare track parameters

    // first track
    dindex sf_idx1 = 0;
    typename bound_track_parameters::vector_type param1;
    getter::element(param1, e_bound_loc0, 0) = 1.;
    getter::element(param1, e_bound_loc1, 0) = 2.;
    getter::element(param1, e_bound_phi, 0) = 0.1;
    getter::element(param1, e_bound_theta, 0) = 0.2;
    getter::element(param1, e_bound_qoverp, 0) = -0.01;
    getter::element(param1, e_bound_time, 0) = 0.1;

    typename bound_track_parameters::covariance_type cov1;

    bound_track_parameters trck1(sf_idx1, param1, cov1);

    // second track
    dindex sf_idx2 = 1;
    typename bound_track_parameters::vector_type param2;
    getter::element(param2, e_bound_loc0, 0) = 4.;
    getter::element(param2, e_bound_loc1, 0) = 20.;
    getter::element(param2, e_bound_phi, 0) = 0.8;
    getter::element(param2, e_bound_theta, 0) = 1.4;
    getter::element(param2, e_bound_qoverp, 0) = 1.;
    getter::element(param2, e_bound_time, 0) = 0.;

    typename bound_track_parameters::covariance_type cov2;

    bound_track_parameters trck2(sf_idx2, param2, cov2);

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

    /// Global position and direction check

    // first track
    auto global_pos1 = trck1.pos(ctx, surfaces, trfs);
    auto u1 = vector::dot(x1, global_pos1 - t1);
    auto v1 = vector::dot(y1, global_pos1 - t1);
    auto w1 = vector::dot(z1, global_pos1 - t1);

    EXPECT_FLOAT_EQ(trck1.local()[0], u1);
    EXPECT_FLOAT_EQ(trck1.local()[1], v1);
    EXPECT_NEAR(0, w1, 1e-6);
    EXPECT_FLOAT_EQ(trck1.dir()[0],
                    std::cos(getter::element(param1, e_bound_phi, 0)) *
                        std::sin(getter::element(param1, e_bound_theta, 0)));
    EXPECT_FLOAT_EQ(trck1.dir()[1],
                    std::sin(getter::element(param1, e_bound_phi, 0)) *
                        std::sin(getter::element(param1, e_bound_theta, 0)));
    EXPECT_FLOAT_EQ(trck1.dir()[2],
                    std::cos(getter::element(param1, e_bound_theta, 0)));

    // second track
    auto global_pos2 = trck2.pos(ctx, surfaces, trfs);
    auto u2 = vector::dot(x2, global_pos2 - t2);
    auto v2 = vector::dot(y2, global_pos2 - t2);
    auto w2 = vector::dot(z2, global_pos2 - t2);

    EXPECT_FLOAT_EQ(trck2.local()[0], u2);
    EXPECT_FLOAT_EQ(trck2.local()[1], v2);
    EXPECT_NEAR(0, w2, 1e-6);
    EXPECT_FLOAT_EQ(trck2.dir()[0],
                    std::cos(getter::element(param2, e_bound_phi, 0)) *
                        std::sin(getter::element(param2, e_bound_theta, 0)));
    EXPECT_FLOAT_EQ(trck2.dir()[1],
                    std::sin(getter::element(param2, e_bound_phi, 0)) *
                        std::sin(getter::element(param2, e_bound_theta, 0)));
    EXPECT_FLOAT_EQ(trck2.dir()[2],
                    std::cos(getter::element(param2, e_bound_theta, 0)));
}

TEST(tools, free_track_parameters) {

    point3 pos = {4., 10., 2.};
    scalar time = 0.1;
    vector3 mom = {10., 20., 30.};
    scalar charge = -1.;

    typename free_track_parameters::vector_type param;
    getter::element(param, e_free_pos0, 0) = pos[0];
    getter::element(param, e_free_pos1, 0) = pos[1];
    getter::element(param, e_free_pos2, 0) = pos[2];
    getter::element(param, e_free_time, 0) = time;
    getter::element(param, e_free_dir0, 0) = mom[0] / getter::norm(mom);
    getter::element(param, e_free_dir1, 0) = mom[1] / getter::norm(mom);
    getter::element(param, e_free_dir2, 0) = mom[2] / getter::norm(mom);
    getter::element(param, e_free_qoverp, 0) = charge / getter::norm(mom);

    typename free_track_parameters::covariance_type cov;

    // first constructor
    free_track_parameters trck1(param, cov);
    EXPECT_FLOAT_EQ(trck1.pos()[0], getter::element(param, e_free_pos0, 0));
    EXPECT_FLOAT_EQ(trck1.pos()[1], getter::element(param, e_free_pos1, 0));
    EXPECT_FLOAT_EQ(trck1.pos()[2], getter::element(param, e_free_pos2, 0));
    EXPECT_FLOAT_EQ(trck1.dir()[0], getter::element(param, e_free_dir0, 0));
    EXPECT_FLOAT_EQ(trck1.dir()[1], getter::element(param, e_free_dir1, 0));
    EXPECT_FLOAT_EQ(trck1.dir()[2], getter::element(param, e_free_dir2, 0));
    EXPECT_FLOAT_EQ(trck1.time(), getter::element(param, e_free_time, 0));
    EXPECT_FLOAT_EQ(trck1.qop(), getter::element(param, e_free_qoverp, 0));
    EXPECT_FLOAT_EQ(trck1.pT(),
                    std::sqrt(std::pow(mom[0], 2) + std::pow(mom[1], 2)));

    // second constructor
    free_track_parameters trck2(pos, time, mom, charge);
    EXPECT_FLOAT_EQ(trck2.pos()[0], pos[0]);
    EXPECT_FLOAT_EQ(trck2.pos()[1], pos[1]);
    EXPECT_FLOAT_EQ(trck2.pos()[2], pos[2]);
    EXPECT_FLOAT_EQ(trck2.dir()[0], mom[0] / getter::norm(mom));
    EXPECT_FLOAT_EQ(trck2.dir()[1], mom[1] / getter::norm(mom));
    EXPECT_FLOAT_EQ(trck2.dir()[2], mom[2] / getter::norm(mom));
    EXPECT_FLOAT_EQ(trck2.time(), time);
    EXPECT_FLOAT_EQ(trck2.qop(), charge / getter::norm(mom));
    EXPECT_FLOAT_EQ(trck2.pT(),
                    std::sqrt(std::pow(mom[0], 2) + std::pow(mom[1], 2)));
}
