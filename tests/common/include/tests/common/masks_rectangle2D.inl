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

constexpr scalar hx{1.f * unit<scalar>::mm};
constexpr scalar hy{9.3f * unit<scalar>::mm};
constexpr scalar hz{0.5f * unit<scalar>::mm};

/// This tests the basic functionality of a rectangle
TEST(mask, rectangle2D) {
    using point_t = typename mask<rectangle2D<>>::loc_point_t;

    point_t p2_in = {0.5f, -9.f};
    point_t p2_edge = {1.f, 9.3f};
    point_t p2_out = {1.5f, -9.f};

    mask<rectangle2D<>> r2{0u, hx, hy};

    ASSERT_NEAR(r2[rectangle2D<>::e_half_x], hx, tol);
    ASSERT_NEAR(r2[rectangle2D<>::e_half_y], hy, tol);

    ASSERT_TRUE(r2.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside(p2_out, 1.f) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = r2.projection_matrix<e_bound_size>();
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j && i < decltype(r2)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = r2.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(hx + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (hx + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], envelope, tol);
}

/// This tests the basic functionality of a cuboid3D
TEST(mask, cuboid3D) {
    using point_t = typename mask<cuboid3D<>>::loc_point_t;

    point_t p2_in = {0.5f, 8.0f, -0.4f};
    point_t p2_edge = {1.f, 9.3f, 0.5f};
    point_t p2_out = {1.5f, -9.f, 0.55f};

    mask<cuboid3D<>> c3{0u, -hx, -hy, -hz, hx, hy, hz};

    ASSERT_NEAR(c3[cuboid3D<>::e_min_x], -hx, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_min_y], -hy, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_min_z], -hz, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_max_x], hx, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_max_y], hy, tol);
    ASSERT_NEAR(c3[cuboid3D<>::e_max_z], hz, tol);

    ASSERT_TRUE(c3.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(c3.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(c3.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c3.is_inside(p2_out, 1.f) == intersection::status::e_inside);

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = c3.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(hx + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -(hz + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (hx + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (hy + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], (hz + envelope), tol);
}
