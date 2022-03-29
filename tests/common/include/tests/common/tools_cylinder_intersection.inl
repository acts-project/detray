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
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/unbound.hpp"
#include "detray/masks/cylinder3.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

// Three-dimensional definitions
using transform3 = __plugin::transform3<detray::scalar>;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();
constexpr scalar not_defined = std::numeric_limits<scalar>::infinity();
constexpr scalar isclose = 1e-5;

// This defines the local frame test suite
TEST(ALGEBRA_PLUGIN, translated_cylinder) {
    // Create a translated cylinder and test untersection
    transform3 shifted(vector3{3., 2., 10.});
    cylinder3<false, cylinder_intersector, unbound, unsigned int>
        cylinder_unbound{4., -10., 10., 0u};
    cylinder_intersector ci;

    // Unbound local frame test
    auto hit_unbound = ci.intersect(shifted, point3{3., 2., 5.},
                                    vector3{1., 0., 0.}, cylinder_unbound);
    ASSERT_TRUE(hit_unbound.status == intersection_status::e_inside);
    ASSERT_NEAR(hit_unbound.p3[0], 7., epsilon);
    ASSERT_NEAR(hit_unbound.p3[1], 2., epsilon);
    ASSERT_NEAR(hit_unbound.p3[2], 5., epsilon);
    ASSERT_TRUE(hit_unbound.p2[0] == not_defined &&
                hit_unbound.p2[1] == not_defined);

    // The same but bound
    cylinder3<false, cylinder_intersector,
              __plugin::cylindrical2<detray::scalar>, unsigned int>
        cylinder_bound{4., -10., 10., 0u};
    auto hit_bound = ci.intersect(shifted, point3{3., 2., 5.},
                                  vector3{1., 0., 0.}, cylinder_bound);
    ASSERT_TRUE(hit_bound.status == intersection_status::e_inside);
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

    // Create a concentric cylinder and test intersection
    scalar r = 4.;
    scalar hz = 10.;
    transform3 identity(vector3{0., 0., 0.});
    cylinder3<false> cylinder{r, -hz, hz, 0u};
    cylinder_intersector ci;
    concentric_cylinder_intersector cci;

    point3 ori = {1., 0.5, 1.};
    point3 dir = vector::normalize(vector3{1., 1., 1.});

    // The same but bound
    auto hit_cylinrical = ci.intersect(identity, ori, dir, cylinder);
    auto hit_cocylindrical = cci.intersect(identity, ori, dir, cylinder);

    ASSERT_TRUE(hit_cylinrical.status == intersection_status::e_inside);
    ASSERT_TRUE(hit_cocylindrical.status == intersection_status::e_inside);
    ASSERT_TRUE(hit_cylinrical.direction == intersection_direction::e_along);
    ASSERT_TRUE(hit_cocylindrical.direction == intersection_direction::e_along);

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
