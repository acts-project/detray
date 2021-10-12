/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "tests/common/read_geometry.hpp"

using namespace detray;

auto [volumes, surfaces, transforms, cylinders, rectangles] = toy_geometry();

// This test check the building of the tml based toy geometry
TEST(ALGEBRA_PLUGIN, toy_geometry) {

    // Check number of geomtery objects
    EXPECT_EQ(volumes.size(), 4);
    EXPECT_EQ(surfaces.size(), 679);
    EXPECT_EQ(transforms.size(), 679);
    EXPECT_EQ(cylinders.size(), 7);
    EXPECT_EQ(rectangles.size(), 672);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
