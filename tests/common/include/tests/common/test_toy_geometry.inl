/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "tests/common/read_geometry.hpp"

using namespace detray;

auto [volumes, surfaces, transforms, cylinders, rectangles] =
    detray_tests::toy_geometry();

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, ray_scan) {

    // Check number of geomtery objects
    EXPECT_EQ(volumes.size(), 1);
    EXPECT_EQ(surfaces.size(), 1);
    EXPECT_EQ(transforms.size(), 1);
    EXPECT_EQ(cylinders.size(), 1);
    EXPECT_EQ(rectangles.size(), 0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
