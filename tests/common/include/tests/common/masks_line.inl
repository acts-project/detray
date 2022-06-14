/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/line.hpp"

using namespace detray;
using point3 = __plugin::point3<scalar>;
using cartesian = __plugin::cartesian2<detray::scalar>;

// This tests the basic function of a line
TEST(mask, line_radial_scope) {

    point3 ln_in{0.1, 0.5, 0};
    point3 ln_edge{1., 50., 0};
    point3 ln_out1{1.2, 0., 0};
    point3 ln_out2{0.1, -51., 0};

    // 50 mm wire with 1 mm radial cell size
    line<> ln{1., 50., 0u};

    ASSERT_FLOAT_EQ(ln[0], 1.);
    ASSERT_FLOAT_EQ(ln[1], 50.);

    ASSERT_TRUE(ln.is_inside<cartesian>(ln_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<cartesian>(ln_edge) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<cartesian>(ln_out1) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside<cartesian>(ln_out2) ==
                intersection::status::e_outside);
}

TEST(mask, line_square_scope) {

    point3 ln_in{1., 0., 0};
    point3 ln_edge{sqrt(2), 0., M_PI / 4};
    point3 ln_out{1.1, 0., 0};

    // 50 mm wire with 1 mm square cell size
    line<ray_line_intersector, cartesian, dindex, true> ln{1., 50., 0u};

    ASSERT_FLOAT_EQ(ln[0], 1.);
    ASSERT_FLOAT_EQ(ln[1], 50.);

    ASSERT_TRUE(ln.is_inside<cartesian>(ln_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<cartesian>(ln_edge, {1e-5, 0}) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<cartesian>(ln_edge, {-1e-5, 0}) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside<cartesian>(ln_out) ==
                intersection::status::e_outside);
}