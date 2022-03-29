/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/ring2.hpp"

using namespace detray;

// This tests the basic function of a rectangle
TEST(mask, ring2) {
    using cartesian = __plugin::cartesian2<detray::scalar>;
    using polar = __plugin::polar2<detray::scalar>;

    point2 p2_pl_in = {0.5, -2.};
    point2 p2_pl_edge = {0., 3.5};
    point2 p2_pl_out = {3.6, 5.};

    point2 p2_c_in = {0.5, -2.};
    point2 p2_c_edge = {0., 3.5};
    point2 p2_c_out = {3.5, 3.5};

    ring2<> r2{0., 3.5, 0u};

    ASSERT_EQ(r2[0], 0.);
    ASSERT_EQ(r2[1], 3.5);

    ASSERT_TRUE(r2.is_inside<polar>(p2_pl_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside<polar>(p2_pl_edge) ==
                intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside<polar>(p2_pl_out) ==
                intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside<polar>(p2_pl_out, 1.2) ==
                intersection::status::e_inside);

    ASSERT_TRUE(r2.is_inside<cartesian>(p2_c_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside<cartesian>(p2_c_edge) ==
                intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside<cartesian>(p2_c_out) ==
                intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside<cartesian>(p2_c_out, 1.45) ==
                intersection::status::e_inside);
}
