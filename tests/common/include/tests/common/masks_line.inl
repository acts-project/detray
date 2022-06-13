/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/line.hpp"

using namespace detray;
using point2 = __plugin::point2<scalar>;

// This tests the basic function of a line
TEST(mask, line) {
    using cartesian = __plugin::cartesian2<detray::scalar>;
    using polar = __plugin::polar2<detray::scalar>;

    point2 ln_p_in = {0.1, 0.0};
    point2 ln_p_edge1 = {1., 2.};
    point2 ln_p_edge2 = {1., 3.};
    point2 ln_p_out = {1.1, 0.};

    point2 ln_c_in = {0.1, 0.5};
    point2 ln_c_edge1 = {1, 0};
    point2 ln_c_edge2 = {0, 1};
    point2 ln_c_out = {0.8, 0.9};

    line<> ln{50., 1., 0u};

    ASSERT_FLOAT_EQ(ln[0], 50.);
    ASSERT_FLOAT_EQ(ln[1], 1.);

    ASSERT_TRUE(ln.is_inside<polar>(ln_p_in) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<polar>(ln_p_edge1) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<polar>(ln_p_edge2) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<polar>(ln_p_out) ==
                intersection::status::e_outside);

    ASSERT_TRUE(ln.is_inside<cartesian>(ln_c_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<cartesian>(ln_c_edge1) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<cartesian>(ln_c_edge2) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside<cartesian>(ln_c_out) ==
                intersection::status::e_outside);
}
