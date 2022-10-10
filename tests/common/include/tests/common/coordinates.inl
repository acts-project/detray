/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/coordinates/coordinates.hpp"
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

// This test cartesian2 coordinate
TEST(test_host_basics, cartesian2) {

    // Preparation work
    const vector3 z = {0., 0., 1.};
    const vector3 x = {1., 0., 0.};
    const point3 t = {2., 3., 4.};
    const transform3 trf(t, z, x);
    const cartesian2<transform3> c2;
    const point3 global1 = {4., 7., 4.};
    const vector3 mom = {1., 2., 3.};
    const vector3 d = vector::normalize(mom);
    const scalar time = 0.1;
    const scalar charge = -1.;
    struct dummy_mask {
    } mask;

    // Global to local transformation
    const point2 local = c2.global_to_local(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], 2., isclose);
    ASSERT_NEAR(local[1], 4., isclose);

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
    ASSERT_NEAR(m.element(bound_vec, 0, 0), 2., isclose);
    ASSERT_NEAR(m.element(bound_vec, 1, 0), 4., isclose);
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
}

// This test cartesian3 coordinate
TEST(test_host_basics, cartesian3) {

    // Preparation work
    const vector3 z = {0., 0., 1.};
    const vector3 x = {1., 0., 0.};
    const point3 t = {2., 3., 4.};
    const transform3 trf(t, z, x);
    const cartesian3<transform3> c3;
    const point3 global1 = {4., 7., 5.};
    const vector3 mom = {1., 2., 3.};
    const vector3 d = vector::normalize(mom);
    struct dummy_mask {
    } mask;

    // Global to local transformation
    const point3 local = c3.global_to_local(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], 2., isclose);
    ASSERT_NEAR(local[1], 4., isclose);
    ASSERT_NEAR(local[2], 1., isclose);

    // Local to global transformation
    const point3 global2 = c3.local_to_global(trf, mask, local, d);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);
}

// This test polar2 coordinate
TEST(test_host_basics, polar2) {

    // Preparation work
    const vector3 z = {0., 0., 1.};
    const vector3 x = {1., 0., 0.};
    const point3 t = {2., 3., 4.};
    const transform3 trf(t, z, x);
    const polar2<transform3> p2;
    const point3 global1 = {4., 7., 4.};
    const vector3 mom = {1., 2., 3.};
    const vector3 d = vector::normalize(mom);
    const scalar time = 0.1;
    const scalar charge = -1.;
    struct dummy_mask {
    } mask;

    // Global to local transformation
    const point2 local = p2.global_to_local(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], std::sqrt(20.), isclose);
    ASSERT_NEAR(local[1], atan2(4., 2.), isclose);

    // Local to global transformation
    const point3 global2 = p2.local_to_global(trf, mask, local, d);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);

    // Free track parameter
    const free_track_parameters<transform3> free_params(global1, time, mom,
                                                        charge);
    const auto free_vec1 = free_params.vector();

    const auto bound_vec = p2.free_to_bound_vector(trf, free_vec1);
    const auto free_vec2 = p2.bound_to_free_vector(trf, mask, bound_vec);

    const matrix_operator m;

    // Check if the bound vector is correct
    ASSERT_NEAR(m.element(bound_vec, 0, 0), std::sqrt(20.), isclose);
    ASSERT_NEAR(m.element(bound_vec, 1, 0), atan2(4., 2.), isclose);
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
}

// This test cylindrical2 coordinate
TEST(test_host_basics, cylindrical2) {

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

    // Define cylinder mask
    struct cylinder_mask {
        scalar r = 0.;
        scalar operator[](dindex) const { return r; }
    };

    const scalar r = 2.;
    const cylinder_mask mask{r};

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
}

// This test cylindrical2 coordinate
TEST(test_host_basics, cylindrical3) {

    // Preparation work
    const vector3 z = {0., 0., 1.};
    const vector3 x = {1., 0., 0.};
    const point3 t = {2., 3., 4.};
    const transform3 trf(t, z, x);
    const cylindrical3<transform3> c3;
    const point3 global1 = {scalar{3.4142136}, scalar{4.4142136}, scalar{9.}};
    const vector3 mom = {1., 2., 3.};
    const vector3 d = vector::normalize(mom);
    struct dummy_mask {
    } mask;

    // Global to local transformation
    const point3 local = c3.global_to_local(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], 2, isclose);
    ASSERT_NEAR(local[1], M_PI_4, isclose);
    ASSERT_NEAR(local[2], 5., isclose);

    // Local to global transformation
    const point3 global2 = c3.local_to_global(trf, mask, local, d);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);
}

// This test line2 coordinate
TEST(test_host_basics, line2) {

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