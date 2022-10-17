/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/coordinates/cylindrical3.hpp"
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

// This test cylindrical3 coordinate
TEST(coordinate, cylindrical3) {

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
