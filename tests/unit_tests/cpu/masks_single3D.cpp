/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/test/types.hpp"

using namespace detray;
using point3_t = test::point3;
using transform3_t = test::transform3;

constexpr scalar tol{1e-7f};

/// This tests the basic functionality of a single value mask (index 0)
GTEST_TEST(detray_masks, single3_0) {
    using point_t = point3_t;

    point_t p3_in = {0.5f, -9.f, 0.f};
    point_t p3_edge = {1.f, 9.3f, 2.f};
    point_t p3_out = {1.5f, -9.8f, 8.f};

    constexpr scalar h0{1.f * unit<scalar>::mm};
    mask<single3D<>> m1_0{0u, -h0, h0};

    ASSERT_NEAR(m1_0[single3D<>::e_lower], -h0, tol);
    ASSERT_NEAR(m1_0[single3D<>::e_upper], h0, tol);

    ASSERT_TRUE(m1_0.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(m1_0.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(m1_0.is_inside(p3_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t0 not t1
    ASSERT_TRUE(m1_0.is_inside(p3_out, 0.6f) == intersection::status::e_inside);

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = m1_0.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(h0 + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (h0 + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], envelope, tol);
}

/// This tests the basic functionality of a single value mask (index 1)
GTEST_TEST(detray_masks, single3_1) {
    using point_t = point3_t;

    point_t p3_in = {0.5f, -9.f, 0.f};
    point_t p3_edge = {1.f, 9.3f, 2.f};
    point_t p3_out = {1.5f, -9.8f, 8.f};

    constexpr scalar h1{9.3f * unit<scalar>::mm};
    mask<single3D<1>> m1_1{0u, -h1, h1};

    ASSERT_NEAR(m1_1[single3D<>::e_lower], -h1, tol);
    ASSERT_NEAR(m1_1[single3D<>::e_upper], h1, tol);

    ASSERT_TRUE(m1_1.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(m1_1.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(m1_1.is_inside(p3_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_1.is_inside(p3_out, 0.6f) == intersection::status::e_inside);

    // Dummy bound track parameter
    bound_track_parameters<transform3_t> bound_params;

    // Check projection matrix
    const auto proj = m1_1.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < decltype(m1_1)::shape::meas_dim; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = m1_1.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(h1 + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (h1 + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], envelope, tol);
}

/// This tests the basic functionality of a single value mask (index 2)
GTEST_TEST(detray_masks, single3_2) {
    using point_t = point3_t;

    point_t p3_in = {0.5f, -9.f, 0.f};
    point_t p3_edge = {1.f, 9.3f, 2.f};
    point_t p3_out = {1.5f, -9.8f, 8.f};

    constexpr scalar h2{2.f * unit<scalar>::mm};
    mask<single3D<2>> m1_2{0u, -h2, h2};

    ASSERT_NEAR(m1_2[single3D<>::e_lower], -h2, tol);
    ASSERT_NEAR(m1_2[single3D<>::e_upper], h2, tol);

    ASSERT_TRUE(m1_2.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(m1_2.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(m1_2.is_inside(p3_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance - take t1 not t1
    ASSERT_TRUE(m1_2.is_inside(p3_out, 6.1f) == intersection::status::e_inside);

    // Dummy bound track parameter
    bound_track_parameters<transform3_t> bound_params;

    // Check projection matrix
    const auto proj = m1_2.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < decltype(m1_2)::shape::meas_dim; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = m1_2.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -(h2 + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], (h2 + envelope), tol);
}
