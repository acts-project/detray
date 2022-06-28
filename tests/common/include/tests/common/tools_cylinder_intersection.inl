/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <cmath>

#include "detray/intersection/concentric_cylinder_intersector.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/masks/cylinder3.hpp"
#include "tests/common/tools/intersectors/helix_cylinder_intersector.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

// Three-dimensional definitions
using transform3 = __plugin::transform3<detray::scalar>;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();
constexpr scalar not_defined = std::numeric_limits<scalar>::infinity();
constexpr scalar isclose = 1e-5;

using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, translated_cylinder) {
    // Create a translated cylinder and test untersection
    const transform3 shifted(vector3{3., 2., 10.});
    cylinder_intersector ci;

    // Test ray
    const point3 ori = {3., 2., 5.};
    const point3 dir = {1., 0., 0.};
    const detail::ray ray(ori, 0., dir, 0.);

    // Check intersection
    cylinder3<cylinder_intersector, __plugin::cylindrical2<detray::scalar>,
              unsigned int>
        cylinder_bound{4., -10., 10., 0u};
    const auto hit_bound = ci(ray, cylinder_bound, shifted)[0];
    ASSERT_TRUE(hit_bound.status == intersection::status::e_inside);
    ASSERT_NEAR(hit_bound.p3[0], 7., epsilon);
    ASSERT_NEAR(hit_bound.p3[1], 2., epsilon);
    ASSERT_NEAR(hit_bound.p3[2], 5., epsilon);
    ASSERT_TRUE(hit_bound.p2[0] != not_defined &&
                hit_bound.p2[1] != not_defined);
    ASSERT_NEAR(hit_bound.p2[0], 0., isclose);
    ASSERT_NEAR(hit_bound.p2[1], -5., isclose);
}

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, concentric_cylinders) {
    // Test ray
    const point3 ori = {1., 0.5, 1.};
    const point3 dir = vector::normalize(vector3{1., 1., 1.});
    const detail::ray ray(ori, 0., dir, 0.);

    // Create a concentric cylinder and test intersection
    const scalar r = 4.;
    const scalar hz = 10.;
    const transform3 identity(vector3{0., 0., 0.});
    cylinder3<> cylinder{r, -hz, hz, 0u};
    cylinder_intersector ci;
    concentric_cylinder_intersector cci;

    // Check intersection
    const auto hit_cylinrical = ci(ray, cylinder, identity)[0];
    const auto hit_cocylindrical = cci(ray, cylinder, identity)[0];

    ASSERT_TRUE(hit_cylinrical.status == intersection::status::e_inside);
    ASSERT_TRUE(hit_cocylindrical.status == intersection::status::e_inside);
    ASSERT_TRUE(hit_cylinrical.direction == intersection::direction::e_along);
    ASSERT_TRUE(hit_cocylindrical.direction ==
                intersection::direction::e_along);

    ASSERT_NEAR(getter::perp(hit_cylinrical.p3), r, isclose);
    ASSERT_NEAR(getter::perp(hit_cocylindrical.p3), r, isclose);

    ASSERT_NEAR(hit_cylinrical.p3[0], hit_cocylindrical.p3[0], isclose);
    ASSERT_NEAR(hit_cylinrical.p3[1], hit_cocylindrical.p3[1], isclose);
    ASSERT_NEAR(hit_cylinrical.p3[2], hit_cocylindrical.p3[2], isclose);
    ASSERT_TRUE(hit_cylinrical.p2[0] != not_defined &&
                hit_cylinrical.p2[1] != not_defined);
    ASSERT_TRUE(hit_cocylindrical.p2[0] != not_defined &&
                hit_cocylindrical.p2[1] != not_defined);
    ASSERT_NEAR(hit_cylinrical.p2[0], hit_cocylindrical.p2[0], isclose);
    ASSERT_NEAR(hit_cylinrical.p2[1], hit_cocylindrical.p2[1], isclose);
}

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, helix_cylinder_intersector) {
    // Create a translated cylinder and test untersection
    const transform3 shifted(vector3{3., 2., 10.});
    helix_cylinder_intersector ci;

    // Test helix
    const point3 pos{3., 2., 5.};
    const vector3 mom{1., 0., 0.};
    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    epsilon * unit_constants::T};
    const detail::helix h({pos, 0, mom, -1}, &B);

    // Check intersection
    cylinder3<cylinder_intersector, __plugin::cylindrical2<detray::scalar>,
              unsigned int>
        cylinder_bound{4., -10., 10., 0u};
    const auto hit_bound = ci(h, cylinder_bound, shifted)[0];
    ASSERT_TRUE(hit_bound.status == intersection::status::e_inside);
    ASSERT_NEAR(hit_bound.p3[0], 7., epsilon);
    ASSERT_NEAR(hit_bound.p3[1], 2., epsilon);
    ASSERT_NEAR(hit_bound.p3[2], 5., epsilon);
    ASSERT_TRUE(hit_bound.p2[0] != not_defined &&
                hit_bound.p2[1] != not_defined);
    ASSERT_NEAR(hit_bound.p2[0], 0., isclose);
    ASSERT_NEAR(hit_bound.p2[1], -5., isclose);
}