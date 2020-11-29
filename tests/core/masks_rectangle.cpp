#include "masks/rectangle2.hpp"
#include "test_defs.hpp"

#include <gtest/gtest.h>

using namespace detray;

// This tests the construction of a surface
TEST(mask, rectangle2)
{
    point2 p2_in = {0.5, -9.};
    point2 p2_edge = {1., 9.3};
    point2 p2_out = {1.5, -9.};

    rectangle2<scalar> r2 = {1., 9.3};

    ASSERT_TRUE(r2(p2_in) == intersection_status::e_inside);
    ASSERT_TRUE(r2(p2_edge) == intersection_status::e_inside);
    ASSERT_TRUE(r2(p2_out) == intersection_status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2(p2_out, 1., 0.5) == intersection_status::e_inside);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
