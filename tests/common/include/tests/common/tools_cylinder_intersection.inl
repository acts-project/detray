/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// Project include(s)
#include "detray/intersection/concentric_cylinder_intersector.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/masks/masks.hpp"
#include "tests/common/tools/intersectors/helix_cylinder_intersector.hpp"

// System include(s)
#include <cmath>
#include <limits>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

namespace {

// Three-dimensional definitions
using transform3_t = __plugin::transform3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;
using ray_t = detray::detail::ray<transform3_t>;
using helix_t = detray::detail::helix<transform3_t>;

constexpr scalar not_defined = std::numeric_limits<scalar>::infinity();
constexpr scalar epsilon = 1e-5f;

constexpr scalar r{4.f};
constexpr scalar hz{10.f};

}  // anonymous namespace

// This checks both solutions of a ray-cylinder intersection
TEST(ALGEBRA_PLUGIN, translated_cylinder) {
    // Create a translated cylinder and test untersection
    const transform3_t shifted(vector3{3.f, 2.f, 10.f});
    cylinder_intersector<transform3_t> ci;

    // Test ray
    const point3 ori = {3.f, 2.f, 5.f};
    const point3 dir = {1.f, 0.f, 0.f};
    const ray_t ray(ori, 0.f, dir, 0.f);

    // Intersect
    mask<cylinder2D<>, unsigned int, transform3_t> cylinder{0u, r, -hz, hz};
    const auto hits_bound = ci(ray, cylinder, shifted);

    // first intersection lies behind the track
    EXPECT_TRUE(hits_bound[0].status == intersection::status::e_inside);
    EXPECT_TRUE(hits_bound[0].direction == intersection::direction::e_opposite);
    EXPECT_NEAR(hits_bound[0].p3[0], -1.f, epsilon);
    EXPECT_NEAR(hits_bound[0].p3[1], 2.f, epsilon);
    EXPECT_NEAR(hits_bound[0].p3[2], 5.f, epsilon);
    ASSERT_TRUE(hits_bound[0].p2[0] != not_defined &&
                hits_bound[0].p2[1] != not_defined);
    // p2[0] = r * phi : 180deg in the opposite direction with r = 4
    EXPECT_NEAR(hits_bound[0].p2[0], 4.f * M_PI, epsilon);
    EXPECT_NEAR(hits_bound[0].p2[1], -5., epsilon);
    EXPECT_NEAR(hits_bound[0].cos_incidence_angle, -1.f, epsilon);

    // second intersection lies in front of the track
    EXPECT_TRUE(hits_bound[1].status == intersection::status::e_inside);
    EXPECT_TRUE(hits_bound[1].direction == intersection::direction::e_along);
    EXPECT_NEAR(hits_bound[1].p3[0], 7.f, epsilon);
    EXPECT_NEAR(hits_bound[1].p3[1], 2.f, epsilon);
    EXPECT_NEAR(hits_bound[1].p3[2], 5.f, epsilon);
    ASSERT_TRUE(hits_bound[1].p2[0] != not_defined &&
                hits_bound[1].p2[1] != not_defined);
    EXPECT_NEAR(hits_bound[1].p2[0], 0.f, epsilon);
    EXPECT_NEAR(hits_bound[1].p2[1], -5.f, epsilon);
    EXPECT_NEAR(hits_bound[1].cos_incidence_angle, 1.f, epsilon);
}

// This checks the inclindence angle calculation for a ray-cylinder intersection
TEST(ALGEBRA_PLUGIN, cylinder_incidence_angle) {
    const transform3_t identity{};
    cylinder_intersector<transform3_t> ci;

    // Test ray
    const point3 ori = {0.f, 1.f, 0.f};
    const point3 dir = {1.f, 0.f, 0.f};
    const ray_t ray(ori, 0.f, dir, 0.f);

    // Intersect
    mask<cylinder2D<>, unsigned int, transform3_t> cylinder{0u, r, -hz, hz};
    const auto hits_bound = ci(ray, cylinder, identity);

    ASSERT_NEAR(hits_bound[0].cos_incidence_angle, -std::sqrt(15.f) / 4.f,
                epsilon);
    ASSERT_NEAR(hits_bound[1].cos_incidence_angle, std::sqrt(15.f) / 4.f,
                epsilon);
}

// This checks the solution of a ray-concentric cylinder intersection against
// those obtained from the portal cylinder intersector.
TEST(ALGEBRA_PLUGIN, concentric_cylinders) {
    // Test ray
    const point3 ori = {1.f, 0.5f, 1.f};
    const point3 dir = vector::normalize(vector3{1.f, 1.f, 1.f});
    const ray_t ray(ori, 0.f, dir, 0.f);

    // Create a concentric cylinder and test intersection
    const transform3_t identity{};
    mask<cylinder2D<>, unsigned int, transform3_t> cylinder{0u, r, -hz, hz};

    cylinder_intersector<transform3_t> ci;
    concentric_cylinder_intersector<transform3_t> cci;

    // Intersect
    const auto hits_cylinrical = ci(ray, cylinder, identity);
    const auto hits_cocylindrical = cci(ray, cylinder, identity);

    ASSERT_TRUE(hits_cylinrical[1].status == intersection::status::e_inside);
    ASSERT_TRUE(hits_cocylindrical[0].status == intersection::status::e_inside);
    ASSERT_TRUE(hits_cylinrical[1].direction ==
                intersection::direction::e_along);
    ASSERT_TRUE(hits_cocylindrical[0].direction ==
                intersection::direction::e_along);

    EXPECT_NEAR(getter::perp(hits_cylinrical[1].p3), r, epsilon);
    EXPECT_NEAR(getter::perp(hits_cocylindrical[0].p3), r, epsilon);

    EXPECT_NEAR(hits_cylinrical[1].p3[0], hits_cocylindrical[0].p3[0], epsilon);
    EXPECT_NEAR(hits_cylinrical[1].p3[1], hits_cocylindrical[0].p3[1], epsilon);
    EXPECT_NEAR(hits_cylinrical[1].p3[2], hits_cocylindrical[0].p3[2], epsilon);
    ASSERT_TRUE(hits_cylinrical[1].p2[0] != not_defined &&
                hits_cylinrical[1].p2[1] != not_defined);
    ASSERT_TRUE(hits_cocylindrical[0].p2[0] != not_defined &&
                hits_cocylindrical[0].p2[1] != not_defined);
    EXPECT_NEAR(hits_cylinrical[1].p2[0], hits_cocylindrical[0].p2[0], epsilon);
    EXPECT_NEAR(hits_cylinrical[1].p2[1], hits_cocylindrical[0].p2[1], epsilon);
}

// This checks the closest solution of a helix-cylinder intersection
TEST(ALGEBRA_PLUGIN, helix_cylinder_intersector) {
    // Create a translated cylinder and test untersection
    const transform3_t shifted(vector3{3.f, 2.f, 10.f});
    detail::helix_cylinder_intersector<transform3_t> hi;

    // Test helix
    const point3 pos{3.f, 2.f, 5.f};
    const vector3 mom{1.f, 0.f, 0.f};
    const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                    epsilon * unit<scalar>::T};
    const helix_t h({pos, 0.f * unit<scalar>::s, mom, -1 * unit<scalar>::e},
                    &B);

    // Intersect
    mask<cylinder2D<false, detail::helix_cylinder_intersector>, unsigned int,
         transform3_t>
        cylinder{0u, r, -hz, hz};
    const auto hits_bound = hi(h, cylinder, shifted);

    // No magnetic field, so the solutions must be the same as for a ray
    ASSERT_TRUE(hits_bound[0].status == intersection::status::e_inside);
    EXPECT_NEAR(hits_bound[0].p3[0], 7.f, epsilon);
    EXPECT_NEAR(hits_bound[0].p3[1], 2.f, epsilon);
    EXPECT_NEAR(hits_bound[0].p3[2], 5.f, epsilon);
    ASSERT_TRUE(hits_bound[0].p2[0] != not_defined &&
                hits_bound[0].p2[1] != not_defined);
    EXPECT_NEAR(hits_bound[0].p2[0], 0.f, epsilon);
    EXPECT_NEAR(hits_bound[0].p2[1], -5.f, epsilon);
    EXPECT_NEAR(hits_bound[0].cos_incidence_angle, 1.f, epsilon);
}