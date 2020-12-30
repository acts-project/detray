/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "grids/axis.hpp"
#include "grids/serializer2.hpp"
#include "utils/indexing.hpp"

#include <gtest/gtest.h>

#include <climits>

using namespace detray;

TEST(grids, serialize_deserialize)
{

    axis::closed<6> r6{-3., 7.};
    axis::circular<12> c12{-3., 3.};

    serializer2 ser2;

    // Serializing
    guaranteed_index test = ser2.serialize<decltype(r6), decltype(c12)>(0u, 0u);
    EXPECT_EQ(test, 0u);
    test = ser2.serialize<decltype(r6), decltype(c12)>(5u, 0u);
    EXPECT_EQ(test, 5u);
    test = ser2.serialize<decltype(r6), decltype(c12)>(0u, 1u);
    EXPECT_EQ(test, 6u);
    test = ser2.serialize<decltype(r6), decltype(c12)>(5u, 2u);
    EXPECT_EQ(test, 17u);

    // Deserialize
    darray<guaranteed_index, 2> expected_array = {0u, 0u};
    darray<guaranteed_index, 2> test_array = ser2.deserialize<decltype(r6), decltype(c12)>(0u);
    EXPECT_EQ(test_array, expected_array);
    expected_array = {5u, 0u};
    test_array = ser2.deserialize<decltype(r6), decltype(c12)>(5u);
    EXPECT_EQ(test_array, expected_array);
    expected_array = {0u, 1u};
    test_array = ser2.deserialize<decltype(r6), decltype(c12)>(6u);
    EXPECT_EQ(test_array, expected_array);
    expected_array = {5u, 2u};
    test_array = ser2.deserialize<decltype(r6), decltype(c12)>(17u);
    EXPECT_EQ(test_array, expected_array);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
