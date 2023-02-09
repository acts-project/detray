/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
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
using transform3_type = __plugin::transform3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;
using scalar_type = typename transform3_type::scalar_type;
using ray_type = detray::detail::ray<transform3_type>;
using helix_type = detray::detail::helix<transform3_type>;

constexpr scalar tol = std::numeric_limits<scalar>::epsilon();
constexpr scalar not_defined = std::numeric_limits<scalar>::infinity();
constexpr scalar isclose{1e-5f};

const scalar r{4.f};
const scalar hz{10.f};

}  // anonymous namespace

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, translated_cylinder) {
    // Create a translated cylinder and test untersection
    const transform3_type shifted(vector3{3.f, 2.f, 10.f});
    cylinder_intersector<transform3_type> ci;

    // Test ray
    const point3 ori = {3.f, 2.f, 5.f};
    const point3 dir = {1.f, 0.f, 0.f};
    const ray_type ray(ori, 0.f, dir, 0.f);

    // Check intersection
    mask<cylinder2D<>, unsigned int, transform3_type> cylinder_bound{0u, r, -hz,
                                                                     hz};
    const auto hit_bound = ci(ray, cylinder_bound, shifted)[0];
    ASSERT_TRUE(hit_bound.status == intersection::status::e_inside);
    ASSERT_NEAR(hit_bound.p3[0], 7.f, tol);
    ASSERT_NEAR(hit_bound.p3[1], 2.f, tol);
    ASSERT_NEAR(hit_bound.p3[2], 5.f, tol);
    ASSERT_TRUE(hit_bound.p2[0] != not_defined &&
                hit_bound.p2[1] != not_defined);
    ASSERT_NEAR(hit_bound.p2[0], 0.f, isclose);
    ASSERT_NEAR(hit_bound.p2[1], -5.f, isclose);
    ASSERT_NEAR(hit_bound.cos_incidence_angle, 1.f, isclose);
}

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, cylinder_incidence_angle) {
    const transform3_type tf(vector3{0.f, 0.f, 0.f});
    cylinder_intersector<transform3_type> ci;

    // Test ray
    const point3 ori = {0.f, 1.f, 0.f};
    const point3 dir = {1.f, 0.f, 0.f};
    const ray_type ray(ori, 0.f, dir, 0.f);

    // Check intersection
    mask<cylinder2D<>, unsigned int, transform3_type> cylinder_bound{0u, r, -hz,
                                                                     hz};
    const auto hit_bound = ci(ray, cylinder_bound, tf)[0];
    ASSERT_NEAR(hit_bound.cos_incidence_angle, 0.25f * std::sqrt(15), isclose);
}

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, concentric_cylinders) {
    // Test ray
    const point3 ori = {1.f, 0.5f, 1.f};
    const point3 dir = vector::normalize(vector3{1.f, 1.f, 1.f});
    const ray_type ray(ori, 0.f, dir, 0.f);

    // Create a concentric cylinder and test intersection
    const transform3_type identity(vector3{0.f, 0.f, 0.f});
    mask<cylinder2D<>, unsigned int, transform3_type> cylinder{0u, r, -hz, hz};
    mask<cylinder2D<false, concentric_cylinder_intersector>, unsigned int,
         transform3_type>
        conc_cylinder{0u, r, -hz, hz};
    cylinder_intersector<transform3_type> ci;
    concentric_cylinder_intersector<transform3_type> cci;

    // Check intersection
    const auto hit_cylinrical = ci(ray, cylinder, identity)[0];
    const auto hit_cocylindrical = cci(ray, conc_cylinder, identity)[0];

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
    const transform3_type shifted(vector3{3.f, 2.f, 10.f});
    helix_cylinder_intersector<transform3_type> hi;

    // Test helix
    const point3 pos{3.f, 2.f, 5.f};
    const vector3 mom{1.f, 0.f, 0.f};
    const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                    tol * unit<scalar>::T};
    const helix_type h({pos, 0.f, mom, -1.f}, &B);

    // Check intersection
    mask<cylinder2D<false, helix_cylinder_intersector>, unsigned int,
         transform3_type>
        cylinder_bound{0u, r, -hz, hz};
    const auto hit_bound = hi(h, cylinder_bound, shifted)[0];
    ASSERT_TRUE(hit_bound.status == intersection::status::e_inside);
    ASSERT_NEAR(hit_bound.p3[0], 7.f, tol);
    ASSERT_NEAR(hit_bound.p3[1], 2.f, tol);
    ASSERT_NEAR(hit_bound.p3[2], 5.f, tol);
    ASSERT_TRUE(hit_bound.p2[0] != not_defined &&
                hit_bound.p2[1] != not_defined);
    ASSERT_NEAR(hit_bound.p2[0], 0.f, isclose);
    ASSERT_NEAR(hit_bound.p2[1], -5.f, isclose);
}