/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/coordinates/cylindrical2.hpp"
#include "detray/masks/masks.hpp"
#include "detray/tracks/tracks.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <climits>

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

// This test cylindrical2 coordinate
TEST(coordinate, cylindrical2) {

    // Preparation work
    const vector3 z = {0., 0., 1.};
    const vector3 x = {1., 0., 0.};
    const point3 t = {2., 3., 4.};
    const transform3 trf(t, z, x);
    const cylindrical2<transform3> c2;
    // Global position on surface
    const point3 global1 = {scalar{3.4142136}, scalar{4.4142136}, scalar{9.}};
    const vector3 mom = {1., 2., 3.};
    const vector3 d = vector::normalize(mom);
    const scalar time = 0.1;
    const scalar charge = -1.;

    const scalar r = 2.;
    const scalar hz = std::numeric_limits<scalar>::infinity();
    mask<cylinder2D<>> mask{0UL, r, -hz, hz};

    // Global to local transformation
    const point2 local = c2.global_to_local(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], r * scalar{M_PI_4}, isclose);
    ASSERT_NEAR(local[1], 5., isclose);

    // Local to global transformation
    const point3 global2 = c2.local_to_global(trf, mask, local, d);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);

    // Free track parameter
    const free_track_parameters<transform3> free_params(global1, time, mom,
                                                        charge);
    const auto free_vec1 = free_params.vector();

    const auto bound_vec = c2.free_to_bound_vector(trf, free_vec1);
    const auto free_vec2 = c2.bound_to_free_vector(trf, mask, bound_vec);

    const matrix_operator m;

    // Check if the bound vector is correct
    ASSERT_NEAR(m.element(bound_vec, 0, 0), r * scalar{M_PI_4}, isclose);
    ASSERT_NEAR(m.element(bound_vec, 1, 0), 5, isclose);
    ASSERT_NEAR(m.element(bound_vec, 2, 0), 1.1071487,
                isclose);  // atan(2)
    ASSERT_NEAR(m.element(bound_vec, 3, 0), 0.64052231,
                isclose);  // atan(sqrt(5)/3)
    ASSERT_NEAR(m.element(bound_vec, 4, 0), -1 / 3.7416574, isclose);
    ASSERT_NEAR(m.element(bound_vec, 5, 0), 0.1, isclose);

    // Check if the same free vector is obtained
    for (int i = 0; i < 8; i++) {
        ASSERT_NEAR(m.element(free_vec1, i, 0), m.element(free_vec2, i, 0),
                    isclose);
    }

    // Normal vector
    const vector3 n =
        c2.normal(trf, mask, free_params.pos(), free_params.dir());
    ASSERT_NEAR(n[0], 1. / std::sqrt(2), isclose);
    ASSERT_NEAR(n[1], 1. / std::sqrt(2), isclose);
    ASSERT_NEAR(n[2], 0., isclose);

    // Test Jacobian transformation
    const matrix_type<6, 6> J =
        c2.free_to_bound_jacobian(trf, mask, free_vec1) *
        c2.bound_to_free_jacobian(trf, mask, bound_vec);

    for (std::size_t i = 0; i < 6; i++) {
        for (std::size_t j = 0; j < 6; j++) {
            if (i == j) {
                EXPECT_NEAR(m.element(J, i, j), 1., isclose);
            } else {
                EXPECT_NEAR(m.element(J, i, j), 0., isclose);
            }
        }
    }
}
