/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/trapezoid2.hpp"

using namespace detray;
using namespace __plugin;

// This tests the basic function of a trapezoid
TEST(mask, trapezoid2) {
    using local_type = __plugin::cartesian2<detray::scalar>;
    using point2 = __plugin::point2<detray::scalar>;

    point2 p2_in = {1., -0.5};
    point2 p2_edge = {2.5, 1.};
    point2 p2_out = {3., 1.5};

    scalar hx_miny = 1.;
    scalar hx_maxy = 3.;
    scalar hy = 2.;
    scalar divisor = 1. / (2. * hy);

    trapezoid2<> t2{hx_miny, hx_maxy, hy, 0u};

    ASSERT_EQ(t2[0], hx_miny);
    ASSERT_EQ(t2[1], hx_maxy);
    ASSERT_EQ(t2[2], hy);
    ASSERT_EQ(t2[3], divisor);

    ASSERT_TRUE(t2.is_inside<local_type>(p2_in) ==
                intersection::status::e_inside);
    ASSERT_TRUE(t2.is_inside<local_type>(p2_edge) ==
                intersection::status::e_inside);
    ASSERT_TRUE(t2.is_inside<local_type>(p2_out) ==
                intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(t2.is_inside<local_type>(p2_out, 1.) ==
                intersection::status::e_inside);
}
