/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/intersection.hpp"
#include "masks/cylinder3.hpp"
#include "masks/single3.hpp"
#include "tools/concentric_cylinder_intersector.hpp"
#include "tools/cylinder_intersector.hpp"
#include "utils/unbound.hpp"

#include <cmath>
#include <climits>

#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

// Three-dimensional definitions
using transform3 = __plugin::transform3;
using vector3 = __plugin::transform3::vector3;
using point3 = __plugin::transform3::point3;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();
constexpr scalar isclose = 1e-5;

// This defines the local frame test suite
TEST(__plugin, translated_cylinder)
{
    // Create a shifted plane
    transform3 shifted(vector3{3., 2., 10.});
    cylinder3<scalar> cylinder = {4., -10., 10.};
    cylinder_intersector ci;

    // Unbound local frame test
    unbound ub;
    auto hit_unbound = ci.intersect(shifted, point3{3., 2., 5.}, vector3{1., 0., 0.}, ub, cylinder);
    ASSERT_TRUE(hit_unbound.status== intersection_status::e_inside);
    ASSERT_NEAR(hit_unbound.point3[0], 7., epsilon);
    ASSERT_NEAR(hit_unbound.point3[1], 2., epsilon);
    ASSERT_NEAR(hit_unbound.point3[2], 5., epsilon);
    ASSERT_TRUE(hit_unbound.point2 == std::nullopt);

    // The same but bound
    __plugin::cylindrical2 cylindrical2;
    auto hit_bound = ci.intersect(shifted, point3{3., 2., 5.}, vector3{1., 0., 0.}, cylindrical2, cylinder);
    ASSERT_TRUE(hit_bound.status== intersection_status::e_inside);
    ASSERT_NEAR(hit_bound.point3[0], 7., epsilon);
    ASSERT_NEAR(hit_bound.point3[1], 2., epsilon);
    ASSERT_NEAR(hit_bound.point3[2], 5., epsilon);
    ASSERT_TRUE(hit_bound.point2 != std::nullopt);
    ASSERT_NEAR(hit_bound.point2.value()[0], 0., isclose);
    ASSERT_NEAR(hit_bound.point2.value()[1], -5., isclose);
}

// This defines the local frame test suite
TEST(__plugin, concentric_cylinders)
{
    using cylinder_surface = surface<transform3, int, int>;

    // Create a shifted plane
    scalar r = 4.;
    scalar hz = 10.;
    transform3 identity(vector3{0., 0., 0.});
    cylinder3<scalar,false> cylinder = {r, -hz, hz};
    cylinder_intersector ci;
    concentric_cylinder_intersector cci;

    point3 ori = {1., 0.5, 1.};
    point3 dir = vector::normalize(vector3{1., 1., 1.});

    // The same but bound
    __plugin::cylindrical2 cylindrical2;
    auto hit_cylinrical = ci.intersect(identity, ori, dir, cylindrical2, cylinder);
    auto hit_cocylinrical = cci.intersect(identity, ori, dir, cylindrical2, cylinder);

    ASSERT_TRUE(hit_cylinrical.status== intersection_status::e_inside);
    ASSERT_TRUE(hit_cocylinrical.status== intersection_status::e_inside);
    ASSERT_TRUE(hit_cylinrical.direction == intersectiondirection::e_along);
    ASSERT_TRUE(hit_cocylinrical.direction == intersectiondirection::e_along);

    ASSERT_NEAR(getter::perp(hit_cylinrical.point3), r, isclose);
    ASSERT_NEAR(getter::perp(hit_cocylinrical.point3), r, isclose);

    ASSERT_NEAR(hit_cylinrical.point3[0], hit_cocylinrical.point3[0], isclose);
    ASSERT_NEAR(hit_cylinrical.point3[1], hit_cocylinrical.point3[1], isclose);
    ASSERT_NEAR(hit_cylinrical.point3[2], hit_cocylinrical.point3[2], isclose);
    ASSERT_TRUE(hit_cylinrical.point2 != std::nullopt);
    ASSERT_TRUE(hit_cocylinrical.point2 != std::nullopt);
    ASSERT_NEAR(hit_cylinrical.point2.value()[0], hit_cocylinrical.point2.value()[0], isclose);
    ASSERT_NEAR(hit_cylinrical.point2.value()[1], hit_cocylinrical.point2.value()[1], isclose);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
