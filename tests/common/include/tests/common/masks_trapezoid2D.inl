/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"

using namespace detray;
using namespace __plugin;

/// This tests the basic functionality of a trapezoid
TEST(mask, trapezoid2D) {
    using point_t = typename mask<trapezoid2D<>>::loc_point_t;

    point_t p2_in = {1., -0.5};
    point_t p2_edge = {2.5, 1.};
    point_t p2_out = {3., 1.5};

    constexpr scalar hx_miny{1. * unit_constants::mm};
    constexpr scalar hx_maxy{3. * unit_constants::mm};
    constexpr scalar hy{2. * unit_constants::mm};
    constexpr scalar divisor{1. / (2. * hy)};

    mask<trapezoid2D<>> t2{0UL, hx_miny, hx_maxy, hy, divisor};

    ASSERT_EQ(t2[trapezoid2D<>::e_half_length_0], hx_miny);
    ASSERT_EQ(t2[trapezoid2D<>::e_half_length_1], hx_maxy);
    ASSERT_EQ(t2[trapezoid2D<>::e_half_length_2], hy);
    ASSERT_EQ(t2[trapezoid2D<>::e_divisor], divisor);

    ASSERT_TRUE(t2.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(t2.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(t2.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(t2.is_inside(p2_out, 1.) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = t2.projection_matrix<e_bound_size>();
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < decltype(t2)::shape::meas_dim; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }

    // Test to_measurement function
    struct test_param {
        using point2 = point_t;
        point_t loc;
        point_t local() const { return loc; }
    } param = {.loc = {1, 2}};

    const auto meas = t2.get_shape().to_measurement(param, {-3, 2});
    ASSERT_EQ(meas, point_t({-2, 4}));
}
