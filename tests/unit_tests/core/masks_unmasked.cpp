/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/unmasked.hpp"
#include "tests/common/test_defs.hpp"

#include <gtest/gtest.h>

using namespace detray;
using namespace plugin;

// This tests the construction of a surface
TEST(mask, unmasked)
{
    point2 p2 = {0.5, -9.};

    unmasked<scalar> u;
    ASSERT_TRUE(u(p2));
    ASSERT_TRUE(u(p2, 1.));
    ASSERT_TRUE(u(p2, 1., 2.));

}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}