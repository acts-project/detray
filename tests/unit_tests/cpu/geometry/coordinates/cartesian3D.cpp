/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/geometry/coordinates/cartesian3D.hpp"

#include "detray/definitions/units.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;
using vector3 = test::vector3;
using transform3 = test::transform3;

const scalar isclose{1e-5f};

// This test cartesian3D coordinate
GTEST_TEST(detray_coordinates, cartesian3D) {

    // Preparation work
    const vector3 z = {0.f, 0.f, 1.f};
    const vector3 x = {1.f, 0.f, 0.f};
    const point3 t = {2.f, 3.f, 4.f};
    const transform3 trf(t, z, x);
    const point3 global1 = {4.f, 7.f, 5.f};
    const vector3 mom = {1.f, 2.f, 3.f};
    const vector3 d = vector::normalize(mom);

    const cartesian3D<test_algebra> c3;

    // Global to local transformation
    const point3 local = c3.global_to_local_3D(trf, global1, d);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], 2.f, isclose);
    ASSERT_NEAR(local[1], 4.f, isclose);
    ASSERT_NEAR(local[2], 1.f, isclose);

    // Local to global transformation
    const point3 global2 = c3.local_to_global(trf, local);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);
}
