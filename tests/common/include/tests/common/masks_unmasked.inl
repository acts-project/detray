/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "masks/unmasked.hpp"

using namespace detray;
using namespace __plugin;

// This tests the construction of a surface
TEST(mask, unmasked) {
    using local_type = __plugin::cartesian2;
    point2 p2 = {0.5, -9.};

    unmasked u;
    ASSERT_TRUE(u.is_inside<local_type>(p2) == e_hit);
    ASSERT_TRUE(u.is_inside<local_type>(p2, true) == e_hit);
    ASSERT_TRUE(u.is_inside<local_type>(p2, false) == e_missed);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
