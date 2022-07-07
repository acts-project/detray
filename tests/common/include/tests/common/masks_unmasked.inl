/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/unmasked.hpp"

using namespace detray;
using point2 = __plugin::point2<scalar>;

// This tests the construction of a surface
TEST(mask, unmasked) {
    point2 p2 = {0.5, -9.};

    unmasked u;
    ASSERT_TRUE(u.is_inside(p2, 0) == intersection::status::e_inside);
}
