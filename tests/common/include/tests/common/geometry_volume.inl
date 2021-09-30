/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <string>

#include "geometry/volume.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(ALGEBRA_PLUGIN, volume) {
    using namespace detray;
    using namespace __plugin;

    using volume = volume<>;
    constexpr bool surface = true;
    constexpr bool portal = false;

    // Check construction, setters and getters
    darray<scalar, 6> bounds = {0., 10., -5., 5., -M_PI, M_PI};
    volume v1 = volume("test_volume", bounds);
    v1.set_index(12345);
    v1.set_surfaces_finder(12);

    ASSERT_TRUE(v1.empty());
    ASSERT_TRUE(v1.index() == 12345);
    ASSERT_TRUE(v1.bounds() == bounds);
    ASSERT_TRUE(v1.surfaces_finder_entry() == 12);

    // Check surface and portal ranges
    dindex_range surface_range{2, 8};
    dindex_range portal_range{20, 24};
    v1.template set_range<surface>(surface_range);
    v1.template set_range<portal>(portal_range);
    ASSERT_TRUE(v1.template range<surface>() == surface_range);
    ASSERT_TRUE(v1.template range<portal>() == portal_range);
    ASSERT_FALSE(v1.empty());
    ASSERT_EQ(v1.template n_objects<surface>(), 6);
    ASSERT_EQ(v1.template n_objects<portal>(), 4);

    // Check copy constructor and eq. operator
    const auto v2 = volume(v1);
    ASSERT_TRUE(v2.index() == 12345);
    ASSERT_TRUE(v2.bounds() == bounds);
    ASSERT_TRUE(v2.surfaces_finder_entry() == 12);
    ASSERT_TRUE(v2.template range<surface>() == surface_range);
    ASSERT_TRUE(v2.template range<portal>() == portal_range);

    ASSERT_TRUE(v2 == v1);
    ASSERT_TRUE(v1 == v2);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
