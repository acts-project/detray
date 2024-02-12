/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/tracks/free_track_parameters.hpp"

#include "detray/test/types.hpp"

// Google Test include(s)
#include <gtest/gtest.h>

using namespace detray;

using vector3 = test::vector3;
using point3 = test::point3;
using transform3 = test::transform3;
using matrix_operator = standard_matrix_operator<scalar>;

constexpr scalar tol{1e-5f};

GTEST_TEST(detray_tracks, free_track_parameters) {

    point3 pos = {4.f, 10.f, 2.f};
    scalar time = 0.1f;
    vector3 mom = {10.f, 20.f, 30.f};
    scalar charge = -1.f;

    typename free_track_parameters<transform3>::vector_type free_vec =
        matrix_operator().template zero<e_free_size, 1>();
    getter::element(free_vec, e_free_pos0, 0u) = pos[0];
    getter::element(free_vec, e_free_pos1, 0u) = pos[1];
    getter::element(free_vec, e_free_pos2, 0u) = pos[2];
    getter::element(free_vec, e_free_time, 0u) = time;
    getter::element(free_vec, e_free_dir0, 0u) = mom[0] / getter::norm(mom);
    getter::element(free_vec, e_free_dir1, 0u) = mom[1] / getter::norm(mom);
    getter::element(free_vec, e_free_dir2, 0u) = mom[2] / getter::norm(mom);
    getter::element(free_vec, e_free_qoverp, 0u) = charge / getter::norm(mom);

    typename free_track_parameters<transform3>::covariance_type free_cov =
        matrix_operator().template zero<e_free_size, e_free_size>();

    // first constructor
    free_track_parameters<transform3> free_param1(free_vec, free_cov);
    EXPECT_NEAR(free_param1.pos()[0],
                getter::element(free_vec, e_free_pos0, 0u), tol);
    EXPECT_NEAR(free_param1.pos()[1],
                getter::element(free_vec, e_free_pos1, 0u), tol);
    EXPECT_NEAR(free_param1.pos()[2],
                getter::element(free_vec, e_free_pos2, 0u), tol);
    EXPECT_NEAR(free_param1.dir()[0],
                getter::element(free_vec, e_free_dir0, 0u), tol);
    EXPECT_NEAR(free_param1.dir()[1],
                getter::element(free_vec, e_free_dir1, 0u), tol);
    EXPECT_NEAR(free_param1.dir()[2],
                getter::element(free_vec, e_free_dir2, 0u), tol);
    EXPECT_NEAR(getter::norm(free_param1.mom()), getter::norm(mom), tol);
    EXPECT_NEAR(free_param1.time(), getter::element(free_vec, e_free_time, 0u),
                tol);
    EXPECT_NEAR(free_param1.qop(), getter::element(free_vec, e_free_qoverp, 0u),
                tol);
    EXPECT_NEAR(free_param1.pT(),
                std::sqrt(std::pow(mom[0], 2.f) + std::pow(mom[1], 2.f)), tol);
    EXPECT_NEAR(free_param1.qopT(), charge / free_param1.pT(), tol);
    EXPECT_NEAR(free_param1.pz(), mom[2], tol);
    EXPECT_NEAR(free_param1.qopz(), charge / free_param1.pz(), tol);
    EXPECT_NEAR(free_param1.mom()[0], free_param1.p() * free_param1.dir()[0],
                tol);
    EXPECT_NEAR(free_param1.mom()[1], free_param1.p() * free_param1.dir()[1],
                tol);
    EXPECT_NEAR(free_param1.mom()[2], free_param1.p() * free_param1.dir()[2],
                tol);

    // second constructor
    free_track_parameters<transform3> free_param2(pos, time, mom, charge);
    EXPECT_NEAR(free_param2.pos()[0], pos[0], tol);
    EXPECT_NEAR(free_param2.pos()[1], pos[1], tol);
    EXPECT_NEAR(free_param2.pos()[2], pos[2], tol);
    EXPECT_NEAR(free_param2.dir()[0], mom[0] / getter::norm(mom), tol);
    EXPECT_NEAR(free_param2.dir()[1], mom[1] / getter::norm(mom), tol);
    EXPECT_NEAR(free_param2.dir()[2], mom[2] / getter::norm(mom), tol);
    EXPECT_NEAR(getter::norm(free_param2.mom()), getter::norm(mom), tol);
    EXPECT_NEAR(free_param2.time(), time, tol);
    EXPECT_NEAR(free_param2.qop(), charge / getter::norm(mom), tol);
    EXPECT_NEAR(free_param2.pT(),
                std::sqrt(std::pow(mom[0], 2.f) + std::pow(mom[1], 2.f)), tol);
    EXPECT_NEAR(free_param2.mom()[0], free_param2.p() * free_param2.dir()[0],
                tol);
    EXPECT_NEAR(free_param2.mom()[1], free_param2.p() * free_param2.dir()[1],
                tol);
    EXPECT_NEAR(free_param2.mom()[2], free_param2.p() * free_param2.dir()[2],
                tol);

    EXPECT_TRUE(free_param2 == free_param1);
}
