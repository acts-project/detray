/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>
#include <climits>

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "grids/axis.hpp"
#include "grids/serializer2.hpp"
#include "utils/indexing.hpp"

using namespace detray;

TEST(grids, serialize_deserialize) {

    axis::regular<> r6{6, -3., 7.};
    axis::circular<> c12{12, -3., 3.};

    serializer2 ser2;

    // Serializing
    dindex test = ser2.serialize(r6, c12, 0u, 0u);
    EXPECT_EQ(test, 0u);
    test = ser2.serialize(r6, c12, 5u, 0u);
    EXPECT_EQ(test, 5u);
    test = ser2.serialize(r6, c12, 0u, 1u);
    EXPECT_EQ(test, 6u);
    test = ser2.serialize(r6, c12, 5u, 2u);
    EXPECT_EQ(test, 17u);

    // Deserialize
    darray<dindex, 2> expected_array = {0u, 0u};
    darray<dindex, 2> test_array = ser2.deserialize(r6, c12, 0u);
    EXPECT_EQ(test_array, expected_array);
    expected_array = {5u, 0u};
    test_array = ser2.deserialize(r6, c12, 5u);
    EXPECT_EQ(test_array, expected_array);
    expected_array = {0u, 1u};
    test_array = ser2.deserialize(r6, c12, 6u);
    EXPECT_EQ(test_array, expected_array);
    expected_array = {5u, 2u};
    test_array = ser2.deserialize(r6, c12, 17u);
    EXPECT_EQ(test_array, expected_array);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
