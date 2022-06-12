/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <cmath>

#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/helix_plane_intersector.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/ray_plane_intersector.hpp"
#include "detray/masks/rectangle2.hpp"
#include "detray/masks/unmasked.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

// Two-dimensional bound frame to surface
__plugin::cartesian2<detray::scalar> cartesian2;

// Three-dimensional definitions
using transform3 = __plugin::transform3<detray::scalar>;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();
constexpr scalar not_defined = std::numeric_limits<scalar>::infinity();
constexpr scalar isclose = 1e-5;

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, translated_plane_ray) {
    // Create a shifted plane
    const transform3 shifted(vector3{3., 2., 10.});
    unmasked<> unmasked_unbound{};
    ray_plane_intersector pi;

    // Test ray
    const point3 pos{2., 1., 0.};
    const vector3 mom{0., 0., 1.};
    const detail::ray r(pos, 0., mom, 0.);

    const auto hit_unbound = pi.intersect(shifted, r, unmasked_unbound);
    ASSERT_TRUE(hit_unbound.status == intersection::status::e_inside);
    ASSERT_TRUE(hit_unbound.direction == intersection::direction::e_along);
    // Global intersection information
    ASSERT_NEAR(hit_unbound.p3[0], 2., epsilon);
    ASSERT_NEAR(hit_unbound.p3[1], 1., epsilon);
    ASSERT_NEAR(hit_unbound.p3[2], 10., epsilon);
    ASSERT_TRUE(hit_unbound.p2[0] == not_defined &&
                hit_unbound.p2[1] == not_defined);

    // The same test but bound to local frame
    unmasked<ray_plane_intersector, __plugin::cartesian2<detray::scalar> >
        unmasked_bound{};
    const auto hit_bound = pi.intersect(shifted, r, unmasked_bound);
    ASSERT_TRUE(hit_bound.status == intersection::status::e_inside);
    // Global intersection information - unchanged
    ASSERT_NEAR(hit_bound.p3[0], 2., epsilon);
    ASSERT_NEAR(hit_bound.p3[1], 1., epsilon);
    ASSERT_NEAR(hit_bound.p3[2], 10., epsilon);
    // Local intersection information
    ASSERT_NEAR(hit_bound.p2[0], -1., epsilon);
    ASSERT_NEAR(hit_bound.p2[1], -1., epsilon);

    // The same test but bound to local frame & masked - inside
    rectangle2<> rect_for_inside{3., 3., 0u};
    const auto hit_bound_inside = pi.intersect(shifted, r, rect_for_inside);
    ASSERT_TRUE(hit_bound_inside.status == intersection::status::e_inside);
    // Global intersection information - unchanged
    ASSERT_NEAR(hit_bound_inside.p3[0], 2., epsilon);
    ASSERT_NEAR(hit_bound_inside.p3[1], 1., epsilon);
    ASSERT_NEAR(hit_bound_inside.p3[2], 10., epsilon);
    // Local intersection infoimation - unchanged
    ASSERT_NEAR(hit_bound_inside.p2[0], -1., epsilon);
    ASSERT_NEAR(hit_bound_inside.p2[1], -1., epsilon);

    // The same test but bound to local frame & masked - outside
    rectangle2<> rect_for_outside{0.5, 3.5, 0u};
    const auto hit_bound_outside = pi.intersect(shifted, r, rect_for_outside);
    ASSERT_TRUE(hit_bound_outside.status == intersection::status::e_outside);
    // Global intersection information - unchanged
    ASSERT_NEAR(hit_bound_outside.p3[0], 2., epsilon);
    ASSERT_NEAR(hit_bound_outside.p3[1], 1., epsilon);
    ASSERT_NEAR(hit_bound_outside.p3[2], 10., epsilon);
    // Local intersection infoimation - unchanged
    ASSERT_NEAR(hit_bound_outside.p2[0], -1., epsilon);
    ASSERT_NEAR(hit_bound_outside.p2[1], -1., epsilon);
}

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, translated_plane_helix) {
    // Create a shifted plane
    const transform3 shifted(vector3{3., 2., 10.});
    unmasked<helix_plane_intersector> unmasked_unbound{};
    helix_plane_intersector pi;

    // Test helix
    const point3 pos{2., 1., 0.};
    const vector3 mom{0., 0., 1.};
    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    epsilon * unit_constants::T};
    const detail::helix h({pos, 0, mom, -1}, &B);

    const auto hit_unbound = pi.intersect(shifted, h, unmasked_unbound);
    ASSERT_TRUE(hit_unbound.status == intersection::status::e_inside);
    ASSERT_TRUE(hit_unbound.direction == intersection::direction::e_along);
    // Global intersection information
    ASSERT_NEAR(hit_unbound.p3[0], 2., epsilon);
    ASSERT_NEAR(hit_unbound.p3[1], 1., epsilon);
    ASSERT_NEAR(hit_unbound.p3[2], 10., epsilon);
    ASSERT_TRUE(hit_unbound.p2[0] == not_defined &&
                hit_unbound.p2[1] == not_defined);

    // The same test but bound to local frame
    unmasked<helix_plane_intersector, __plugin::cartesian2<detray::scalar> >
        unmasked_bound{};
    const auto hit_bound = pi.intersect(shifted, h, unmasked_bound);
    ASSERT_TRUE(hit_bound.status == intersection::status::e_inside);
    // Global intersection information - unchanged
    ASSERT_NEAR(hit_bound.p3[0], 2., epsilon);
    ASSERT_NEAR(hit_bound.p3[1], 1., epsilon);
    ASSERT_NEAR(hit_bound.p3[2], 10., epsilon);
    // Local intersection information
    ASSERT_NEAR(hit_bound.p2[0], -1., epsilon);
    ASSERT_NEAR(hit_bound.p2[1], -1., epsilon);

    // The same test but bound to local frame & masked - inside
    rectangle2<helix_plane_intersector> rect_for_inside{3., 3., 0u};
    const auto hit_bound_inside = pi.intersect(shifted, h, rect_for_inside);
    ASSERT_TRUE(hit_bound_inside.status == intersection::status::e_inside);
    // Global intersection information - unchanged
    ASSERT_NEAR(hit_bound_inside.p3[0], 2., epsilon);
    ASSERT_NEAR(hit_bound_inside.p3[1], 1., epsilon);
    ASSERT_NEAR(hit_bound_inside.p3[2], 10., epsilon);
    // Local intersection infoimation - unchanged
    ASSERT_NEAR(hit_bound_inside.p2[0], -1., epsilon);
    ASSERT_NEAR(hit_bound_inside.p2[1], -1., epsilon);

    // The same test but bound to local frame & masked - outside
    rectangle2<helix_plane_intersector> rect_for_outside{0.5, 3.5, 0u};
    const auto hit_bound_outside = pi.intersect(shifted, h, rect_for_outside);
    ASSERT_TRUE(hit_bound_outside.status == intersection::status::e_outside);
    // Global intersection information - unchanged
    ASSERT_NEAR(hit_bound_outside.p3[0], 2., epsilon);
    ASSERT_NEAR(hit_bound_outside.p3[1], 1., epsilon);
    ASSERT_NEAR(hit_bound_outside.p3[2], 10., epsilon);
    // Local intersection infoimation - unchanged
    ASSERT_NEAR(hit_bound_outside.p2[0], -1., epsilon);
    ASSERT_NEAR(hit_bound_outside.p2[1], -1., epsilon);
}
