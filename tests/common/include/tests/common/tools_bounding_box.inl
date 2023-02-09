/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/tools/bounding_box.hpp"

using namespace detray;
using namespace __plugin;

namespace {

// cuboid
constexpr scalar hx{1.f * unit<scalar>::mm};
constexpr scalar hy{9.3f * unit<scalar>::mm};
constexpr scalar hz{0.5f * unit<scalar>::mm};

// cylinder
constexpr scalar r{3.f * unit<scalar>::mm};

// envelope around wrapped object (scalor in percent)
constexpr scalar envelope{1.01f};

// test tolerance
constexpr scalar tol{1e-7f};

}  // anonymous namespace

/// This tests the basic functionality cuboid axis aligned bounding box
TEST(tools, cuboid_bounding_box) {
    using point_t = typename mask<cylinder3D>::loc_point_t;

    point_t p2_in = {0.5f, 8.0f, -0.4f};
    point_t p2_edge = {1.f, 9.3f, 0.5f};
    point_t p2_out = {1.5f, -9.f, 0.55f};

    mask<cuboid3D<>> c3{0u, -hx, hx, -hy, hy, -hz, hz};

    axis_aligned_bounding_box<> aabb{c3, 0u};

    // Id of this instance
    ASSERT_EQ(aabb.id(), 0u);

    // Test the bounds
    const auto bounds = aabb.bounds();
    ASSERT_NEAR(bounds[cuboid3D<>::e_min_x], -hx * envelope, tol);
    ASSERT_NEAR(bounds[cuboid3D<>::e_min_y], -hy * envelope, tol);
    ASSERT_NEAR(bounds[cuboid3D<>::e_min_z], -hz * envelope, tol);
    ASSERT_NEAR(bounds[cuboid3D<>::e_max_x], hx * envelope, tol);
    ASSERT_NEAR(bounds[cuboid3D<>::e_max_y], hy * envelope, tol);
    ASSERT_NEAR(bounds[cuboid3D<>::e_max_z], hz * envelope, tol);

    ASSERT_TRUE(aabb.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(aabb.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(aabb.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(aabb.is_inside(p2_out, 1.) == intersection::status::e_inside);
}

/// This tests the basic functionality cylindrical axis aligned bounding box
TEST(tools, cyl_bounding_box) {
    using point_t = typename mask<cylinder3D>::loc_point_t;

    point_t p3_in = {r, 0.f, -1.f};
    point_t p3_edge = {0.f, r, hz};
    point_t p3_out = {static_cast<scalar>(r / std::sqrt(2.)),
                      static_cast<scalar>(r / std::sqrt(2.)), 4.5};
    point_t p3_off = {1.f, 1.f, -9.f};

    mask<cylinder3D> c3{
        0u,  0.f, r, static_cast<scalar>(-M_PI), static_cast<scalar>(M_PI),
        -hz, hz};

    axis_aligned_bounding_box<cylinder3D> aabb{c3, 0u};

    // Id of this instance
    ASSERT_EQ(aabb.id(), 0u);

    // Test the bounds
    const auto bounds = aabb.bounds();
    ASSERT_NEAR(bounds[cylinder3D::e_min_r], 0.f, tol);
    ASSERT_NEAR(bounds[cylinder3D::e_max_r], r * envelope, tol);
    ASSERT_NEAR(bounds[cylinder3D::e_min_phi], -M_PI * envelope, tol);
    ASSERT_NEAR(bounds[cylinder3D::e_max_phi], M_PI * envelope, tol);
    ASSERT_NEAR(bounds[cylinder3D::e_min_z], -hz * envelope, tol);
    ASSERT_NEAR(bounds[cylinder3D::e_max_z], hz * envelope, tol);

    ASSERT_TRUE(aabb.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(aabb.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(aabb.is_inside(p3_out) == intersection::status::e_outside);
    ASSERT_TRUE(aabb.is_inside(p3_off) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(aabb.is_inside(p3_out, 0.6) == intersection::status::e_inside);
}