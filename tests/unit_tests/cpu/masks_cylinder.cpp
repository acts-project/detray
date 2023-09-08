/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"
#include "detray/test/types.hpp"

using namespace detray;
using point3_t = test::point3;
using transform3_t = test::transform3;

constexpr scalar tol{1e-7f};

constexpr scalar r{3.f * unit<scalar>::mm};
constexpr scalar hz{4.f * unit<scalar>::mm};

/// This tests the basic functionality of a 2D cylinder
GTEST_TEST(detray_masks, cylinder2D) {
    using point_t = point3_t;

    point_t p2_in = {r, -1.f, r};
    point_t p2_edge = {r, hz, r};
    point_t p2_out = {3.5f, 4.5f, 3.5f};

    mask<cylinder2D<>> c{0u, r, -hz, hz};

    ASSERT_NEAR(c[cylinder2D<>::e_r], r, tol);
    ASSERT_NEAR(c[cylinder2D<>::e_n_half_z], -hz, tol);
    ASSERT_NEAR(c[cylinder2D<>::e_p_half_z], hz, tol);

    ASSERT_TRUE(c.is_inside(p2_in));
    ASSERT_TRUE(c.is_inside(p2_edge));
    ASSERT_FALSE(c.is_inside(p2_out));
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c.is_inside(p2_out, 0.6f));

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = c.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -(hz + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], (hz + envelope), tol);
}

/// This tests the basic functionality of a 3D cylinder
GTEST_TEST(detray_masks, cylinder3D) {
    using point_t = point3_t;

    point_t p3_in = {r, 0.f, -1.f};
    point_t p3_edge = {0.f, r, hz};
    point_t p3_out = {r * constant<scalar>::inv_sqrt2,
                      r * constant<scalar>::inv_sqrt2, 4.5f};
    point_t p3_off = {1.f, 1.f, -9.f};

    mask<cylinder3D> c{
        0u, 0.f, -constant<scalar>::pi, -hz, r, constant<scalar>::pi, hz};

    ASSERT_NEAR(c[cylinder3D::e_min_r], 0.f, tol);
    ASSERT_NEAR(c[cylinder3D::e_max_r], r, tol);
    ASSERT_NEAR(c[cylinder3D::e_min_phi], -constant<scalar>::pi, tol);
    ASSERT_NEAR(c[cylinder3D::e_max_phi], constant<scalar>::pi, tol);
    ASSERT_NEAR(c[cylinder3D::e_min_z], -hz, tol);
    ASSERT_NEAR(c[cylinder3D::e_max_z], hz, tol);

    ASSERT_TRUE(c.is_inside(p3_in));
    ASSERT_TRUE(c.is_inside(p3_edge));
    ASSERT_FALSE(c.is_inside(p3_out));
    ASSERT_FALSE(c.is_inside(p3_off));
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c.is_inside(p3_out, 0.6f));

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = c.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -(hz + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], (hz + envelope), tol);
}
