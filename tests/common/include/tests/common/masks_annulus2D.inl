/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/masks.hpp"

using namespace detray;
using namespace __plugin;

/// This tests the basic functionality of a stereo annulus
TEST(mask, annulus2D) {
    using point_t = typename mask<annulus2D<>>::loc_point_t;

    scalar minR = 7.2;
    scalar maxR = 12.0;
    scalar minPhi = 0.74195;
    scalar maxPhi = 1.33970;
    point_t offset = {-2., 2.};

    // points in cartesian module frame
    point_t p2_in = {7., 7.};
    point_t p2_out1 = {5., 5.};
    point_t p2_out2 = {10., 3.};
    point_t p2_out3 = {10., 10.};
    point_t p2_out4 = {4., 10.};

    auto toStripFrame = [&](const point_t& xy) -> point_t {
        auto shifted = xy + offset;
        scalar r = getter::perp(shifted);
        scalar phi = getter::phi(shifted);
        return point_t{r, phi};
    };

    mask<annulus2D<>> ann2{0UL,    minR,      maxR,      minPhi,
                           maxPhi, offset[0], offset[1], 0.f};

    ASSERT_FLOAT_EQ(ann2[0], static_cast<scalar>(7.2));
    ASSERT_FLOAT_EQ(ann2[1], static_cast<scalar>(12.0));
    ASSERT_FLOAT_EQ(ann2[2], static_cast<scalar>(0.74195));
    ASSERT_FLOAT_EQ(ann2[3], static_cast<scalar>(1.33970));
    ASSERT_FLOAT_EQ(ann2[4], static_cast<scalar>(-2.0));
    ASSERT_FLOAT_EQ(ann2[5], static_cast<scalar>(2.0));
    ASSERT_FLOAT_EQ(ann2[6], static_cast<scalar>(0.));

    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_in)) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out1)) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out2)) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out3)) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out4)) ==
                intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out1), 1.3) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out4), 0.07) ==
                intersection::status::e_inside);
}
