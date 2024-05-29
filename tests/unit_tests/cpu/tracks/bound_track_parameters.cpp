/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/tracks/bound_track_parameters.hpp"

#include "detray/test/common/types.hpp"

// Google Test include(s)
#include <gtest/gtest.h>

using namespace detray;

using algebra_t = test::algebra;
using vector3 = test::vector3;
using point3 = test::point3;
using matrix_operator = test::matrix_operator;

constexpr scalar tol{1e-5f};

GTEST_TEST(detray_tracks, bound_track_parameters) {

    /// Declare track parameters

    // first track
    dindex sf_idx1 = 0u;
    typename bound_track_parameters<algebra_t>::vector_type bound_vec1 =
        matrix_operator().template zero<e_bound_size, 1>();
    getter::element(bound_vec1, e_bound_loc0, 0u) = 1.f;
    getter::element(bound_vec1, e_bound_loc1, 0u) = 2.f;
    getter::element(bound_vec1, e_bound_phi, 0u) = 0.1f;
    getter::element(bound_vec1, e_bound_theta, 0u) = 0.2f;
    getter::element(bound_vec1, e_bound_qoverp, 0u) = -0.01f;
    getter::element(bound_vec1, e_bound_time, 0u) = 0.1f;

    typename bound_track_parameters<algebra_t>::covariance_type bound_cov1 =
        matrix_operator().template zero<e_bound_size, e_bound_size>();

    bound_track_parameters<algebra_t> bound_param1(
        geometry::barcode{}.set_index(sf_idx1), bound_vec1, bound_cov1);
    EXPECT_NEAR(bound_param1.pT(),
                1.f /
                    std::abs(getter::element(bound_vec1, e_bound_qoverp, 0u)) *
                    std::sin(getter::element(bound_vec1, e_bound_theta, 0u)),
                tol);
    EXPECT_NEAR(bound_param1.qopT(), -1.f / bound_param1.pT(), tol);
    EXPECT_NEAR(bound_param1.pz(),
                1.f /
                    std::abs(getter::element(bound_vec1, e_bound_qoverp, 0u)) *
                    std::cos(getter::element(bound_vec1, e_bound_theta, 0u)),
                tol);
    EXPECT_NEAR(bound_param1.qopz(), -1.f / bound_param1.pz(), tol);

    // second track
    dindex sf_idx2 = 1u;
    typename bound_track_parameters<algebra_t>::vector_type bound_vec2 =
        matrix_operator().template zero<e_bound_size, 1>();
    getter::element(bound_vec2, e_bound_loc0, 0u) = 4.f;
    getter::element(bound_vec2, e_bound_loc1, 0u) = 20.f;
    getter::element(bound_vec2, e_bound_phi, 0u) = 0.8f;
    getter::element(bound_vec2, e_bound_theta, 0u) = 1.4f;
    getter::element(bound_vec2, e_bound_qoverp, 0u) = 1.f;
    getter::element(bound_vec2, e_bound_time, 0u) = 0.f;

    typename bound_track_parameters<algebra_t>::covariance_type bound_cov2 =
        matrix_operator().template zero<e_bound_size, e_bound_size>();

    bound_track_parameters<algebra_t> bound_param2(
        geometry::barcode{}.set_index(sf_idx2), bound_vec2, bound_cov2);
    bound_track_parameters<algebra_t> bound_param3(
        geometry::barcode{}.set_index(sf_idx2), bound_vec2, bound_cov2);

    /// Check the elements

    // first track
    EXPECT_NEAR(bound_param1.bound_local()[0],
                getter::element(bound_vec1, e_bound_loc0, 0u), tol);
    EXPECT_NEAR(bound_param1.bound_local()[1],
                getter::element(bound_vec1, e_bound_loc1, 0u), tol);
    EXPECT_NEAR(bound_param1.phi(),
                getter::element(bound_vec1, e_bound_phi, 0u), tol);
    EXPECT_NEAR(bound_param1.theta(),
                getter::element(bound_vec1, e_bound_theta, 0u), tol);
    EXPECT_NEAR(bound_param1.qop(),
                getter::element(bound_vec1, e_bound_qoverp, 0u), tol);
    EXPECT_NEAR(bound_param1.charge(), -1.f, tol);
    EXPECT_NEAR(bound_param1.time(),
                getter::element(bound_vec1, e_bound_time, 0u), tol);
    EXPECT_NEAR(bound_param1.mom()[0],
                bound_param1.p() * std::sin(bound_param1.theta()) *
                    std::cos(bound_param1.phi()),
                tol);
    EXPECT_NEAR(bound_param1.mom()[1],
                bound_param1.p() * std::sin(bound_param1.theta()) *
                    std::sin(bound_param1.phi()),
                tol);
    EXPECT_NEAR(bound_param1.mom()[2],
                bound_param1.p() * std::cos(bound_param1.theta()), tol);

    // second track
    EXPECT_NEAR(bound_param2.bound_local()[0],
                getter::element(bound_vec2, e_bound_loc0, 0u), tol);
    EXPECT_NEAR(bound_param2.bound_local()[1],
                getter::element(bound_vec2, e_bound_loc1, 0u), tol);
    EXPECT_NEAR(bound_param2.phi(),
                getter::element(bound_vec2, e_bound_phi, 0u), tol);
    EXPECT_NEAR(bound_param2.theta(),
                getter::element(bound_vec2, e_bound_theta, 0u), tol);
    EXPECT_NEAR(bound_param2.qop(),
                getter::element(bound_vec2, e_bound_qoverp, 0u), tol);
    EXPECT_NEAR(bound_param2.charge(), 1.f, tol);
    EXPECT_NEAR(bound_param2.time(),
                getter::element(bound_vec2, e_bound_time, 0u), tol);
    EXPECT_NEAR(bound_param2.mom()[0],
                bound_param2.p() * std::sin(bound_param2.theta()) *
                    std::cos(bound_param2.phi()),
                tol);
    EXPECT_NEAR(bound_param2.mom()[1],
                bound_param2.p() * std::sin(bound_param2.theta()) *
                    std::sin(bound_param2.phi()),
                tol);
    EXPECT_NEAR(bound_param2.mom()[2],
                bound_param2.p() * std::cos(bound_param2.theta()), tol);

    EXPECT_TRUE(!(bound_param2 == bound_param1));
    EXPECT_TRUE(bound_param2 == bound_param3);
}
