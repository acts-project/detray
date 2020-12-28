/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "grids/axis.hpp"
#include "utils/indexing.hpp"

#include <gtest/gtest.h>

#include <climits>

using namespace detray;

TEST(grids, regular_closed_axis)
{

    axis::closed<10> ten_bins{-3., 7.};
    // N bins
    EXPECT_EQ(ten_bins.bins, 10u);
    // Axis bin access
    EXPECT_EQ(ten_bins.bin(-4.), 0u);
    EXPECT_EQ(ten_bins.bin(2.), 5u);
    EXPECT_EQ(ten_bins.bin(8.), 9u);
    // Axis range access
    guaranteed_range expected_range = {4u, 6u};
    EXPECT_EQ(ten_bins.range(2., 1), expected_range);
    expected_range = {0u, 8u};
    EXPECT_EQ(ten_bins.range(1., 4), expected_range);
    expected_range = {3u, 9u};
    EXPECT_EQ(ten_bins.range(5., 5), expected_range);
    // Axis sequence access
    guaranteed_sequence expected_zone = {4u, 5u, 6u};
    EXPECT_EQ(ten_bins.zone(2., 1), expected_zone);
    expected_zone = {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u};
    EXPECT_EQ(ten_bins.zone(1., 4), expected_zone);

}

TEST(grids, regular_circular_axis)
{
    scalar epsilon = 10 * std::numeric_limits<scalar>::epsilon();

    // Let's say 36 modules, but with 4 directly at 0, pi/2, pi, -pi2
    scalar half_module = 2 * M_PI_2 / 72;
    axis::circular<36> full_pi = {-M_PI + half_module, M_PI - half_module};
    // N bins
    EXPECT_EQ(full_pi.bins, 36u);
    // Axis bin access
    EXPECT_EQ(full_pi.bin(M_PI - epsilon), 0u);
    EXPECT_EQ(full_pi.bin(M_PI + epsilon), 0u);
    EXPECT_EQ(full_pi.bin(0), 18u);
    // Axis range access
    guaranteed_range expected = {4u, 6u};
    //EXPECT_EQ(ten_bins.zone(2., 1), expected);
    //expected = {0u, 8u};
    //EXPECT_EQ(ten_bins.zone(1., 4), expected);
    //expected = { 3u, 9u};
    //EXPECT_EQ(ten_bins.zone(5., 5), expected);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
