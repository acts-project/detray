/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/coordinates/line2.hpp"
#include "detray/tracks/tracks.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;
using point2 = __plugin::point2<scalar>;
using point3 = __plugin::point3<scalar>;
using vector3 = __plugin::vector3<scalar>;
using transform3 = __plugin::transform3<detray::scalar>;
using matrix_operator = typename transform3::matrix_actor;
using size_type = typename matrix_operator::size_ty;
template <size_type ROWS, size_type COLS>
using matrix_type = typename matrix_operator::template matrix_type<ROWS, COLS>;

const scalar isclose = 1e-5;

TEST(coordinate, line2_case1) {

    // Preparation work
    vector3 z = {1., 1., 1.};
    z = vector::normalize(z);
    vector3 x = {1., 0., -1.};
    x = vector::normalize(x);
    const point3 t = {0., 0., 0.};
    const transform3 trf(t, z, x);
    const line2<transform3> l2;
    const point3 global1 = {1, 1.5, 0.5};
    const vector3 mom = {0., 1., 1.};
    const vector3 d = vector::normalize(mom);
    const scalar time = 0.1;
    const scalar charge = -1.;
    struct dummy_mask {
    } mask;

    // Global to local transformation
    const point2 local = l2.global_to_local(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], -1. / sqrt(2), isclose);
    ASSERT_NEAR(local[1], sqrt(3), isclose);

    // Local to global transformation
    const point3 global2 = l2.local_to_global(trf, mask, local, d);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);

    // Free track parameter
    const free_track_parameters<transform3> free_params(global1, time, mom,
                                                        charge);
    const auto free_vec1 = free_params.vector();

    const auto bound_vec = l2.free_to_bound_vector(trf, free_vec1);
    const auto free_vec2 = l2.bound_to_free_vector(trf, mask, bound_vec);

    const matrix_operator m;

    // Check if the bound vector is correct
    ASSERT_NEAR(m.element(bound_vec, 0, 0), -1. / sqrt(2), isclose);
    ASSERT_NEAR(m.element(bound_vec, 1, 0), sqrt(3), isclose);
    ASSERT_NEAR(m.element(bound_vec, 2, 0), M_PI_2, isclose);
    ASSERT_NEAR(m.element(bound_vec, 3, 0), M_PI_4, isclose);
    ASSERT_NEAR(m.element(bound_vec, 4, 0), -1. / sqrt(2), isclose);
    ASSERT_NEAR(m.element(bound_vec, 5, 0), 0.1, isclose);

    // Check if the same free vector is obtained
    for (int i = 0; i < 8; i++) {
        ASSERT_NEAR(m.element(free_vec1, i, 0), m.element(free_vec2, i, 0),
                    isclose);
    }

    // Test Jacobian transformation
    const matrix_type<6, 6> J =
        l2.free_to_bound_jacobian(trf, mask, free_vec1) *
        l2.bound_to_free_jacobian(trf, mask, bound_vec);

    for (std::size_t i = 0; i < 6; i++) {
        for (std::size_t j = 0; j < 6; j++) {
            if (i == j) {
                EXPECT_NEAR(m.element(J, i, j), 1., isclose);
            } else {
                EXPECT_NEAR(m.element(J, i, j), 0., isclose);
            }
        }
    }

    const auto free_to_bound = l2.free_to_bound_jacobian(trf, mask, free_vec1);
    const auto bound_to_free = l2.bound_to_free_jacobian(trf, mask, bound_vec);

    for (std::size_t i = 0; i < 8; i++) {
        for (std::size_t j = 0; j < 8; j++) {
            printf("%f ", m.element(J, i, j));
        }
        printf("\n");
    }
    printf("\n");

    for (std::size_t i = 0; i < 6; i++) {
        for (std::size_t j = 0; j < 8; j++) {
            printf("%f ", m.element(free_to_bound, i, j));
        }
        printf("\n");
    }
    printf("\n");

    for (std::size_t i = 0; i < 8; i++) {
        for (std::size_t j = 0; j < 6; j++) {
            printf("%f ", m.element(bound_to_free, i, j));
        }
        printf("\n");
    }
}

TEST(coordinate, line2_case2) {

    // Preparation work
    vector3 z = {1., 2., 3.};
    z = vector::normalize(z);
    vector3 x = {2., -4., 2.};
    x = vector::normalize(x);
    const point3 t = {0., 0., 0.};
    const transform3 trf(t, z, x);
    const line2<transform3> l2;
    const point2 local1 = {1, 2};
    const vector3 mom = {1., 6., -2.};
    const vector3 d = vector::normalize(mom);
    const scalar time = 0.1;
    const scalar charge = -1.;
    struct dummy_mask {
    } mask;

    // local to global transformation
    const point3 global = l2.local_to_global(trf, mask, local1, d);

    // global to local transform
    const point2 local2 = l2.global_to_local(trf, global, d);

    // Check if the same local position is obtained
    ASSERT_NEAR(local1[0], local2[0], isclose);
    ASSERT_NEAR(local1[1], local2[1], isclose);

    // Free track parameter
    const free_track_parameters<transform3> free_params(global, time, mom,
                                                        charge);
    const auto free_vec = free_params.vector();
    const auto bound_vec = l2.free_to_bound_vector(trf, free_vec);

    // Test Jacobian transformation
    const matrix_type<6, 6> J = l2.free_to_bound_jacobian(trf, mask, free_vec) *
                                l2.bound_to_free_jacobian(trf, mask, bound_vec);

    const matrix_operator m;

    for (std::size_t i = 0; i < 6; i++) {
        for (std::size_t j = 0; j < 6; j++) {
            if (i == j) {
                EXPECT_NEAR(m.element(J, i, j), 1., isclose);
            } else {
                EXPECT_NEAR(m.element(J, i, j), 0., isclose);
            }
        }
    }

    const auto free_to_bound = l2.free_to_bound_jacobian(trf, mask, free_vec);
    const auto bound_to_free = l2.bound_to_free_jacobian(trf, mask, bound_vec);

    for (std::size_t i = 0; i < 8; i++) {
        for (std::size_t j = 0; j < 8; j++) {
            printf("%f ", m.element(J, i, j));
        }
        printf("\n");
    }
    printf("\n");

    for (std::size_t i = 0; i < 6; i++) {
        for (std::size_t j = 0; j < 8; j++) {
            printf("%f ", m.element(free_to_bound, i, j));
        }
        printf("\n");
    }
    printf("\n");

    for (std::size_t i = 0; i < 8; i++) {
        for (std::size_t j = 0; j < 6; j++) {
            printf("%f ", m.element(bound_to_free, i, j));
        }
        printf("\n");
    }
}

// This test line2 coordinate
TEST(coordinate, line2_case3) {
    /*
    // Preparation work
    const vector3 z = {0., 0., 1.};
    const vector3 x = {1., 0., 0.};
    const point3 t = {2., 3., 4.};
    const transform3 trf(t, z, x);
    const line2<transform3> l2;
    const point3 global1 = {3., 3., 9.};
    const vector3 mom = {0., 2., 0.};
    const vector3 d = vector::normalize(mom);
    const scalar time = 0.1;
    const scalar charge = -1.;
    struct dummy_mask {
    } mask;

    // Global to local transformation
    const point2 local = l2.global_to_local(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], -1., isclose);
    ASSERT_NEAR(local[1], 5., isclose);

    // Local to global transformation
    const point3 global2 = l2.local_to_global(trf, mask, local, d);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);

    // Free track parameter
    const free_track_parameters<transform3> free_params(global1, time, mom,
                                                        charge);
    const auto free_vec1 = free_params.vector();

    const auto bound_vec = l2.free_to_bound_vector(trf, free_vec1);
    const auto free_vec2 = l2.bound_to_free_vector(trf, mask, bound_vec);

    const matrix_operator m;

    // Check if the bound vector is correct
    ASSERT_NEAR(m.element(bound_vec, 0, 0), -1., isclose);
    ASSERT_NEAR(m.element(bound_vec, 1, 0), 5, isclose);
    ASSERT_NEAR(m.element(bound_vec, 2, 0), M_PI_2, isclose);
    ASSERT_NEAR(m.element(bound_vec, 3, 0), M_PI_2, isclose);
    ASSERT_NEAR(m.element(bound_vec, 4, 0), -0.5, isclose);
    ASSERT_NEAR(m.element(bound_vec, 5, 0), 0.1, isclose);

    // Check if the same free vector is obtained
    for (int i = 0; i < 8; i++) {
        ASSERT_NEAR(m.element(free_vec1, i, 0), m.element(free_vec2, i, 0),
                    isclose);
    }

    // Test Jacobian transformation
    const matrix_type<6, 6> J =
        l2.free_to_bound_jacobian(trf, mask, free_vec1) *
        l2.bound_to_free_jacobian(trf, mask, bound_vec);

    for (std::size_t i = 0; i < 6; i++) {
        for (std::size_t j = 0; j < 6; j++) {
            if (i == j) {
                EXPECT_NEAR(m.element(J, i, j), 1., isclose);
            } else {
                EXPECT_NEAR(m.element(J, i, j), 0., isclose);
            }
        }
    }
    */
}

// This test line2 coordinate
TEST(coordinate, line2_case4) {
    /*
    // Preparation work
    vector3 z = {0., 1., 1.};
    z = vector::normalize(z);
    const vector3 x = {1., 0., 0.};
    // const point3 t = {0., 0., 0.};
    const point3 t = {0., -1., -1.};
    const transform3 trf(t, z, x);
    const line2<transform3> l2;
    const point3 global1 = {1., 0., 0.};
    const vector3 mom = {0., 2., 0.};
    const vector3 d = vector::normalize(mom);
    const scalar time = 0.1;
    const scalar charge = -1.;
    struct dummy_mask {
    } mask;

    // Global to local transformation
    const point2 local = l2.global_to_local(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], -1., isclose);
    // ASSERT_NEAR(local[1], 0., isclose);
    ASSERT_NEAR(local[1], std::sqrt(2), isclose);

    // Local to global transformation
    const point3 global2 = l2.local_to_global(trf, mask, local, d);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);

    // Free track parameter
    const free_track_parameters<transform3> free_params(global1, time, mom,
                                                        charge);
    const auto free_vec1 = free_params.vector();

    const auto bound_vec = l2.free_to_bound_vector(trf, free_vec1);
    const auto free_vec2 = l2.bound_to_free_vector(trf, mask, bound_vec);

    const matrix_operator m;

    // Check if the bound vector is correct
    ASSERT_NEAR(m.element(bound_vec, 0, 0), -1., isclose);
    // ASSERT_NEAR(m.element(bound_vec, 1, 0), 0, isclose);
    ASSERT_NEAR(m.element(bound_vec, 1, 0), std::sqrt(2), isclose);
    ASSERT_NEAR(m.element(bound_vec, 2, 0), M_PI_2, isclose);
    ASSERT_NEAR(m.element(bound_vec, 3, 0), M_PI_2, isclose);
    ASSERT_NEAR(m.element(bound_vec, 4, 0), -0.5, isclose);
    ASSERT_NEAR(m.element(bound_vec, 5, 0), 0.1, isclose);

    // Check if the same free vector is obtained
    for (int i = 0; i < 8; i++) {
        ASSERT_NEAR(m.element(free_vec1, i, 0), m.element(free_vec2, i, 0),
                    isclose);
    }

    // Test Jacobian transformation
    const matrix_type<6, 6> J =
        l2.free_to_bound_jacobian(trf, mask, free_vec1) *
        l2.bound_to_free_jacobian(trf, mask, bound_vec);

    for (std::size_t i = 0; i < 6; i++) {
        for (std::size_t j = 0; j < 6; j++) {
            if (i == j) {
                EXPECT_NEAR(m.element(J, i, j), 1., isclose);
            } else {
                EXPECT_NEAR(m.element(J, i, j), 0., isclose);
            }
        }
    }

    const auto free_to_bound = l2.free_to_bound_jacobian(trf, mask, free_vec1);
    const auto bound_to_free = l2.bound_to_free_jacobian(trf, mask, bound_vec);

    for (std::size_t i = 0; i < 8; i++) {
        for (std::size_t j = 0; j < 8; j++) {
            printf("%f ", m.element(J, i, j));
        }
        printf("\n");
    }
    printf("\n");

    for (std::size_t i = 0; i < 6; i++) {
        for (std::size_t j = 0; j < 8; j++) {
            printf("%f ", m.element(free_to_bound, i, j));
        }
        printf("\n");
    }
    printf("\n");

    for (std::size_t i = 0; i < 8; i++) {
        for (std::size_t j = 0; j < 6; j++) {
            printf("%f ", m.element(bound_to_free, i, j));
        }
        printf("\n");
    }
    */
}