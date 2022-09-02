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

const scalar isclose = 1e-5;

// This test line2 coordinate
TEST(coordinate, line2) {

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
}