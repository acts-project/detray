/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/annulus.hpp"

#include <gtest/gtest.h>

using namespace detray;
using namespace __plugin;

// This tests the basic function of a rectangle
TEST(mask, annulus)
{
    //using scalar = __plugin::scalar;
    using point2_c    = __plugin::cartesian2::point2;
    using point2_pol  = __plugin::polar2::point2;

    scalar minR = 7.2;
    scalar maxR = 12.0;
    scalar minPhi = 0.74195;
    scalar maxPhi = 1.33970;
    point2_c offset = {-2., 2.};

    // points in cartesian module frame
    point2_c p2_in   = {7., 7.};
    point2_c p2_out1 = {5., 5.};
    point2_c p2_out2 = {10., 3.};
    point2_c p2_out3 = {10., 10.};
    point2_c p2_out4 = {4., 10.};

    auto toStripFrame = [&](const point2_pol& xy) -> point2_pol {
        auto shifted = xy + offset;
        scalar r   = getter::perp(shifted);
        scalar phi = getter::phi(shifted);
        return point2_pol{r, phi};
    };

    annulus<scalar> ann2 = {minR, maxR, minPhi, maxPhi, offset[0], offset[1]};

    ASSERT_EQ(ann2[0], 7.2);
    ASSERT_EQ(ann2[1], 12.0);
    ASSERT_EQ(ann2[2], 0.74195);
    ASSERT_EQ(ann2[3], 1.33970);
    ASSERT_EQ(ann2[4], -2.0);
    ASSERT_EQ(ann2[5], 2.0);
    ASSERT_EQ(ann2[6], 0.);

    ASSERT_TRUE(ann2(toStripFrame(p2_in), 0., 0.) == intersection_status::e_inside);
    ASSERT_TRUE(ann2(toStripFrame(p2_out1), 0., 0.) == intersection_status::e_outside);
    ASSERT_TRUE(ann2(toStripFrame(p2_out2), 0., 0.) == intersection_status::e_outside);
    ASSERT_TRUE(ann2(toStripFrame(p2_out3), 0., 0.) == intersection_status::e_outside);
    ASSERT_TRUE(ann2(toStripFrame(p2_out4), 0., 0.) == intersection_status::e_outside);
    // Move outside point inside using a tolerance
    //ASSERT_TRUE(ann2(toStripFrame(p2_out1), t0, t1, t2, t3) == intersection_status::e_outside);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
