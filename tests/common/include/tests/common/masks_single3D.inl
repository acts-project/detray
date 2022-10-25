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

/// This tests the basic functionality of a single value mask (index 0)
TEST(mask, single3_0) {
    using point_t = typename mask<single3D<>>::loc_point_t;

    point_t p3_in = {0.5, -9., 0.};
    point_t p3_edge = {1., 9.3, 2.};
    point_t p3_out = {1.5, -9.8, 8.};

    constexpr scalar h0{1. * unit_constants::mm};
    mask<single3D<>> m1_0{0UL, -h0, h0};

    ASSERT_FLOAT_EQ(m1_0[single3D<>::e_lower], -h0);
    ASSERT_FLOAT_EQ(m1_0[single3D<>::e_upper], h0);

    ASSERT_TRUE(m1_0.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(m1_0.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(m1_0.is_inside(p3_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t0 not t1
    ASSERT_TRUE(m1_0.is_inside(p3_out, 0.6) == intersection::status::e_inside);
}

/// This tests the basic functionality of a single value mask (index 1)
TEST(mask, single3_1) {
    using point_t = typename mask<single3D<1>>::loc_point_t;

    point_t p3_in = {0.5, -9., 0.};
    point_t p3_edge = {1., 9.3, 2.};
    point_t p3_out = {1.5, -9.8, 8.};

    constexpr scalar h1{9.3 * unit_constants::mm};
    mask<single3D<1>> m1_1{0UL, -h1, h1};

    ASSERT_FLOAT_EQ(m1_1[single3D<>::e_lower], -h1);
    ASSERT_FLOAT_EQ(m1_1[single3D<>::e_upper], h1);

    ASSERT_TRUE(m1_1.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(m1_1.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(m1_1.is_inside(p3_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_1.is_inside(p3_out, 0.6) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = m1_1.projection_matrix<e_bound_size>();
    for (std::size_t i = 0; i < decltype(m1_1)::shape::meas_dim; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }
}

/// This tests the basic functionality of a single value mask (index 2)
TEST(mask, single3_2) {
    using point_t = typename mask<single3D<2>>::loc_point_t;

    point_t p3_in = {0.5, -9., 0.};
    point_t p3_edge = {1., 9.3, 2.};
    point_t p3_out = {1.5, -9.8, 8.};

    constexpr scalar h2{2. * unit_constants::mm};
    mask<single3D<2>> m1_2{0UL, -h2, h2};

    ASSERT_FLOAT_EQ(m1_2[single3D<>::e_lower], -h2);
    ASSERT_FLOAT_EQ(m1_2[single3D<>::e_upper], h2);

    ASSERT_TRUE(m1_2.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(m1_2.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(m1_2.is_inside(p3_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_2.is_inside(p3_out, 6.1) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = m1_2.projection_matrix<e_bound_size>();
    for (std::size_t i = 0; i < decltype(m1_2)::shape::meas_dim; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }
}
