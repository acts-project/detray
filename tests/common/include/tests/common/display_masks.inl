/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <cmath>

#include "detray/masks/masks.hpp"
#include "style/styles.hpp"
#include "view/draw.hpp"
#include "view/views.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

using transform3 = __plugin::transform3;
using vector3 = __plugin::vector3;
using point3 = __plugin::point3;

TEST(display, annulus2) {
    detray::global_xy_view gxy;

    // First rectangle
    transform3 transform{};
    annulus2<> annulus = {7.2, 12.0, 0.74195, 1.33970, 0., -2., 2.};

    color c = {0.2, 0.8, 0.6, 0.9};
    display(false);
    draw_mask(annulus, transform, style{c}, gxy);
    save("annulus.png");
}

TEST(display, rectangle2) {
    detray::global_xy_view gxy;

    // First rectangle
    transform3 transform{};
    rectangle2<> rectangle = {3., 4.};

    color c = {0.5, 0.2, 0.6, 0.9};
    display(false);
    draw_mask(rectangle, transform, style{c}, gxy);
    save("rectangle.png");
}

TEST(display, trapezoid2) {
    detray::global_xy_view gxy;

    // First rectangle
    transform3 transform{};
    trapezoid2<> trapezoid = {3., 4., 5.};

    color c = {0.5, 0.4, 0.4, 0.9};
    display(false);
    draw_mask(trapezoid, transform, style{c}, gxy);
    save("trapezoid.png");
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
