/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

/// TODO: These headers don't exist??

// #include <gtest/gtest.h>

// #include <climits>
// #include <cmath>

// #include "detray/masks/masks.hpp"
// #include "detray/style/styles.hpp"
// #include "detray/view/draw.hpp"
// #include "detray/view/views.hpp"
// #include "detray/test/types.hpp"

// /// @note test has to be defined with a preprocessor command
// using namespace detray;

// using transform3 = test::transform3;
// using vector3 = test::vector3;
// using point3 = test::point3;

// GTEST_TEST(detray_core, annulus2) {
//     detray::global_xy_view gxy;

//     // First rectangle
//     transform3 transform{};
//     annulus2<> annulus = {7.2, 12.0, 0.74195, 1.33970, 0., -2., 2.};

//     color c = {0.2, 0.8, 0.6, 0.9};
//     display(false);
//     draw_mask(annulus, transform, style{c}, gxy);
//     save("annulus.png");
// }

// GTEST_TEST(detray_core, rectangle2) {
//     detray::global_xy_view gxy;

//     // First rectangle
//     transform3 transform{};
//     rectangle2<> rectangle = {3., 4.};

//     color c = {0.5, 0.2, 0.6, 0.9};
//     display(false);
//     draw_mask(rectangle, transform, style{c}, gxy);
//     save("rectangle.png");
// }

// GTEST_TEST(detray_core, trapezoid2) {
//     detray::global_xy_view gxy;

//     // First rectangle
//     transform3 transform{};
//     trapezoid2<> trapezoid = {3., 4., 5.};

//     color c = {0.5, 0.4, 0.4, 0.9};
//     display(false);
//     draw_mask(trapezoid, transform, style{c}, gxy);
//     save("trapezoid.png");
// }
