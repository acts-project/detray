/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"

using namespace detray;
using namespace __plugin;

constexpr scalar tol{1e-7f};

/// This tests the basic functionality of a trapezoid
TEST(mask, trapezoid2D) {
    using point_t = typename mask<trapezoid2D<>>::loc_point_t;

    point_t p2_in = {1.f, -0.5f};
    point_t p2_edge = {2.5f, 1.f};
    point_t p2_out = {3.f, 1.5f};

    constexpr scalar hx_miny{1.f * unit<scalar>::mm};
    constexpr scalar hx_maxy{3.f * unit<scalar>::mm};
    constexpr scalar hy{2.f * unit<scalar>::mm};
    constexpr scalar divisor{1.f / (2.f * hy)};

    mask<trapezoid2D<>> t2{0u, hx_miny, hx_maxy, hy, divisor};

    ASSERT_NEAR(t2[trapezoid2D<>::e_half_length_0], hx_miny, tol);
    ASSERT_NEAR(t2[trapezoid2D<>::e_half_length_1], hx_maxy, tol);
    ASSERT_NEAR(t2[trapezoid2D<>::e_half_length_2], hy, tol);
    ASSERT_NEAR(t2[trapezoid2D<>::e_divisor], divisor, tol);

    ASSERT_TRUE(t2.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(t2.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(t2.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(t2.is_inside(p2_out, 1.) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = t2.projection_matrix<e_bound_size>();
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j && i < decltype(t2)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = t2.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(hx_maxy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (hx_maxy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], envelope, tol);
}
