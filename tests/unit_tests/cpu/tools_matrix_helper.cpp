/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/utils/matrix_helper.hpp"
#include "detray/test/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;
using transform3 = test::transform3;
using matrix_operator = standard_matrix_operator<scalar>;
using vector3 = typename transform3::vector3;

constexpr scalar tolerance = 1e-6f;

GTEST_TEST(detray_core, column_wise_cross) {
    auto P = matrix_operator().template zero<3, 3>();

    getter::element(P, 0u, 0u) = 0.f;
    getter::element(P, 0u, 1u) = 1.f;
    getter::element(P, 0u, 2u) = 2.f;
    getter::element(P, 1u, 0u) = 3.f;
    getter::element(P, 1u, 1u) = 4.f;
    getter::element(P, 1u, 2u) = 5.f;
    getter::element(P, 2u, 0u) = 6.f;
    getter::element(P, 2u, 1u) = 7.f;
    getter::element(P, 2u, 2u) = 8.f;

    const vector3 u{1.f, 2.f, 3.f};

    const auto Q = matrix_helper<matrix_operator>().column_wise_cross(P, u);

    EXPECT_NEAR(getter::element(Q, 0u, 0u), -3.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 1u, 0u), 6.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 2u, 0u), -3.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 0u, 1u), -2.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 1u, 1u), 4.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 2u, 1u), -2.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 0u, 2u), -1.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 1u, 2u), 2.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 2u, 2u), -1.f, tolerance);
}

GTEST_TEST(detray_core, column_wise_multiply) {

    auto P = matrix_operator().template zero<3, 3>();

    getter::element(P, 0u, 0u) = 0.f;
    getter::element(P, 0u, 1u) = 1.f;
    getter::element(P, 0u, 2u) = 2.f;
    getter::element(P, 1u, 0u) = 3.f;
    getter::element(P, 1u, 1u) = 4.f;
    getter::element(P, 1u, 2u) = 5.f;
    getter::element(P, 2u, 0u) = 6.f;
    getter::element(P, 2u, 1u) = 7.f;
    getter::element(P, 2u, 2u) = 8.f;

    const vector3 u{1.f, 2.f, 3.f};

    const auto Q = matrix_helper<matrix_operator>().column_wise_multiply(P, u);

    EXPECT_NEAR(getter::element(Q, 0u, 0u), 0.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 0u, 1u), 1.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 0u, 2u), 2.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 1u, 0u), 6.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 1u, 1u), 8.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 1u, 2u), 10.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 2u, 0u), 18.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 2u, 1u), 21.f, tolerance);
    EXPECT_NEAR(getter::element(Q, 2u, 2u), 24.f, tolerance);
}

GTEST_TEST(detray_core, cross_matrix) {

    const vector3 u{1.f, 2.f, 3.f};
    const vector3 v{3.f, 4.f, 5.f};

    const auto u_cross = matrix_helper<matrix_operator>().cross_matrix(u);
    const auto v_cross = matrix_helper<matrix_operator>().cross_matrix(v);

    EXPECT_NEAR(getter::element(u_cross, 0u, 0u), 0.f, tolerance);
    EXPECT_NEAR(getter::element(u_cross, 0u, 1u), -3.f, tolerance);
    EXPECT_NEAR(getter::element(u_cross, 0u, 2u), 2.f, tolerance);
    EXPECT_NEAR(getter::element(u_cross, 1u, 0u), 3.f, tolerance);
    EXPECT_NEAR(getter::element(u_cross, 1u, 1u), 0.f, tolerance);
    EXPECT_NEAR(getter::element(u_cross, 1u, 2u), -1.f, tolerance);
    EXPECT_NEAR(getter::element(u_cross, 2u, 0u), -2.f, tolerance);
    EXPECT_NEAR(getter::element(u_cross, 2u, 1u), 1.f, tolerance);
    EXPECT_NEAR(getter::element(u_cross, 2u, 2u), 0.f, tolerance);

    // [u]_cross * v = [v]_cross^T * u
    const vector3 u_cross_v = u_cross * v;
    const vector3 v_cross_u = matrix_operator().transpose(v_cross) * u;

    EXPECT_NEAR(u_cross_v[0], -2.f, tolerance);
    EXPECT_NEAR(u_cross_v[1], 4.f, tolerance);
    EXPECT_NEAR(u_cross_v[2], -2.f, tolerance);
    EXPECT_NEAR(u_cross_v[0], v_cross_u[0], tolerance);
    EXPECT_NEAR(u_cross_v[1], v_cross_u[1], tolerance);
    EXPECT_NEAR(u_cross_v[2], v_cross_u[2], tolerance);
}

GTEST_TEST(detray_core, outer_product) {

    const vector3 u{1.f, 2.f, 3.f};
    const vector3 v{3.f, 4.f, 5.f};

    const auto m33 = matrix_helper<matrix_operator>().outer_product(u, v);

    EXPECT_NEAR(getter::element(m33, 0u, 0u), 3.f, tolerance);
    EXPECT_NEAR(getter::element(m33, 0u, 1u), 4.f, tolerance);
    EXPECT_NEAR(getter::element(m33, 0u, 2u), 5.f, tolerance);
    EXPECT_NEAR(getter::element(m33, 1u, 0u), 6.f, tolerance);
    EXPECT_NEAR(getter::element(m33, 1u, 1u), 8.f, tolerance);
    EXPECT_NEAR(getter::element(m33, 1u, 2u), 10.f, tolerance);
    EXPECT_NEAR(getter::element(m33, 2u, 0u), 9.f, tolerance);
    EXPECT_NEAR(getter::element(m33, 2u, 1u), 12.f, tolerance);
    EXPECT_NEAR(getter::element(m33, 2u, 2u), 15.f, tolerance);
}