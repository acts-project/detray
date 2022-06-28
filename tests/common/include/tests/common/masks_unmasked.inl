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
    using local_type = __plugin::cartesian2<detray::scalar>;
    point2 p2 = {0.5, -9.};

    unmasked u;
    ASSERT_TRUE(u.is_inside<local_type>(p2, 0) ==
                intersection::status::e_inside);
    /*
    ASSERT_TRUE(u.is_inside<local_type>() ==
                intersection::status::e_inside);
    ASSERT_TRUE(u.is_inside<local_type>(p2, false) ==
                intersection::status::e_outside);
    */
}
