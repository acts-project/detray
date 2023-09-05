/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/coordinates/line2D.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using point2 = test::point2;
using point3 = test::point3;
using vector3 = test::vector3;
using transform3 = test::transform3;
using matrix_operator = typename transform3::matrix_actor;
using size_type = typename matrix_operator::size_ty;
template <size_type ROWS, size_type COLS>
using matrix_type = typename matrix_operator::template matrix_type<ROWS, COLS>;

constexpr scalar isclose{1e-5f};

GTEST_TEST(detray_coordinates, line2_case1) {

    // Preparation work
    vector3 z = {1.f, 1.f, 1.f};
    z = vector::normalize(z);
    vector3 x = {1.f, 0.f, -1.f};
    x = vector::normalize(x);
    const point3 t = {0.f, 0.f, 0.f};
    const transform3 trf(t, z, x);
    const line2D<ALGEBRA_PLUGIN<test::scalar>> l2;
    const point3 global1 = {1.f, 1.5f, 0.5f};
    const vector3 mom = {0.f, 1.f, 1.f};
    const vector3 d = vector::normalize(mom);
    const scalar time{0.1f};
    const scalar charge{-1.f};
    struct dummy_mask {
    } mask;

    // Global to local transformation
    const point3 local = l2.global_to_local(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], -constant<scalar>::inv_sqrt2, isclose);
    ASSERT_NEAR(local[1], std::sqrt(3.f), isclose);

    // Local to global transformation
    const point3 global2 = l2.local_to_global(trf, local);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);

    // Normal vector
    const vector3 n = l2.normal(trf);
    ASSERT_EQ(n, z);

    // Free track parameter
    /*const free_track_parameters<transform3> free_params(global1, time, mom,
                                                        charge);
    const auto free_vec1 = free_params.vector();

    const auto bound_vec = l2.free_to_bound_vector(trf, free_vec1);
    const auto free_vec2 = l2.bound_to_free_vector(trf, mask, bound_vec);

    const matrix_operator m;

    // Check if the bound vector is correct
    ASSERT_NEAR(m.element(bound_vec, 0u, 0u), -constant<scalar>::inv_sqrt2,
                isclose);
    ASSERT_NEAR(m.element(bound_vec, 1u, 0u), std::sqrt(3.f), isclose);
    ASSERT_NEAR(m.element(bound_vec, 2u, 0u), constant<scalar>::pi_2, isclose);
    ASSERT_NEAR(m.element(bound_vec, 3u, 0u), constant<scalar>::pi_4, isclose);
    ASSERT_NEAR(m.element(bound_vec, 4u, 0u), -constant<scalar>::inv_sqrt2,
                isclose);
    ASSERT_NEAR(m.element(bound_vec, 5u, 0u), 0.1f, isclose);

    // Check if the same free vector is obtained
    for (unsigned int i = 0u; i < 8u; i++) {
        ASSERT_NEAR(m.element(free_vec1, i, 0u), m.element(free_vec2, i, 0u),
                    isclose);
    }

    // Test Jacobian transformation
    const matrix_type<6, 6> J = l2.free_to_bound_jacobian(trf, free_vec1) *
                                l2.bound_to_free_jacobian(trf, mask, bound_vec);

    for (unsigned int i = 0u; i < 6u; i++) {
        for (unsigned int j = 0u; j < 6u; j++) {
            if (i == j) {
                EXPECT_NEAR(m.element(J, i, j), 1.f, isclose);
            } else {
                EXPECT_NEAR(m.element(J, i, j), 0.f, isclose);
            }
        }
    }*/
}

GTEST_TEST(detray_coordinates, line2_case2) {

    // Preparation work
    vector3 z = {1.f, 2.f, 3.f};
    z = vector::normalize(z);
    vector3 x = {2.f, -4.f, 2.f};
    x = vector::normalize(x);
    const point3 t = {0.f, 0.f, 0.f};
    const transform3 trf(t, z, x);
    const line2D<ALGEBRA_PLUGIN<test::scalar>> l2;
    const point2 local1 = {1.f, 2.f};
    const vector3 mom = {1.f, 6.f, -2.f};
    const vector3 d = vector::normalize(mom);
    const scalar time{0.1f};
    const scalar charge{-1.f};
    struct dummy_mask {
    } mask;

    // local to global transformation
    const point3 global = l2.bound_local_to_global(trf, mask, local1, d);

    // global to local transform
    const point3 local2 = l2.global_to_local(trf, global, d);

    // Check if the same local position is obtained
    ASSERT_NEAR(local1[0], local2[0], isclose);
    ASSERT_NEAR(local1[1], local2[1], isclose);

    // Free track parameter
    /*const free_track_parameters<transform3> free_params(global, time, mom,
                                                        charge);
    const auto free_vec = free_params.vector();
    const auto bound_vec = l2.free_to_bound_vector(trf, free_vec);

    // Test Jacobian transformation
    const matrix_type<6, 6> J = l2.free_to_bound_jacobian(trf, free_vec) *
                                l2.bound_to_free_jacobian(trf, mask, bound_vec);

    const matrix_operator m;

    for (unsigned int i = 0u; i < 6u; i++) {
        for (unsigned int j = 0u; j < 6u; j++) {
            if (i == j) {
                EXPECT_NEAR(m.element(J, i, j), 1.f, isclose);
            } else {
                EXPECT_NEAR(m.element(J, i, j), 0.f, isclose);
            }
        }
    }*/
}
