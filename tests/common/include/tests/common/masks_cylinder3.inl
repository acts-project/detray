/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/cylinder3.hpp"

using namespace detray;
using namespace __plugin;

// This tests the basic function of a rectangle
TEST(mask, cylinder3) {
    using local_type = __plugin::transform3<detray::scalar>;
    using point3 = __plugin::point3<detray::scalar>;

    scalar r = 3.;
    scalar hz = 4.;

    point3 p3_in = {r, 0., -1.};
    point3 p3_edge = {0., r, hz};
    point3 p3_out = {static_cast<scalar>(r / std::sqrt(2.)),
                     static_cast<scalar>(r / std::sqrt(2.)), 4.5};
    point3 p3_off = {1., 1., -9.};

    // Test radius to be on surface, too
    cylinder3<cylinder_intersector, local_type, dindex, true> c{r, -hz, hz, 0u};

    ASSERT_EQ(c[0], r);
    ASSERT_EQ(c[1], -hz);
    ASSERT_EQ(c[2], hz);

    ASSERT_TRUE(c.is_inside<local_type>(p3_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(c.is_inside<local_type>(p3_edge) ==
                intersection::status::e_inside);
    ASSERT_TRUE(c.is_inside<local_type>(p3_out) ==
                intersection::status::e_outside);
    ASSERT_TRUE(c.is_inside<local_type>(p3_off) ==
                intersection::status::e_missed);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c.is_inside<local_type>(p3_out, {0., 0.6}) ==
                intersection::status::e_inside);
}
