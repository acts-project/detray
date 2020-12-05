/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/surface.hpp"
#include "core/intersection.hpp"
#include "masks/cylinder3.hpp"
#include "masks/single3.hpp"
#include "tools/concentric_cylinder_intersector.hpp"
#include "tools/cylinder_intersector.hpp"
#include "utils/unbound.hpp"

#include <cmath>
#include <climits>

#include <gtest/gtest.h>

/// @note plugin has to be defined with a preprocessor command
using namespace detray;

// Three-dimensional definitions
using transform3 = plugin::transform3;
using vector3 = plugin::transform3::vector3;
using point3 = plugin::transform3::point3;
using context = plugin::transform3::context;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();
constexpr scalar isclose = 1e-5;

// This defines the local frame test suite
TEST(plugin, translated_cylinder)
{
    context ctx;
    using cylinder_surface = surface<transform3, int, int>;

    // Create a shifted plane
    transform3 shifted(vector3{3., 2., 10.}, ctx);
    cylinder_surface shifted_cylinder(std::move(shifted), 1, 1);
    cylinder3<scalar> cylinder = {4., 10.};
    cylinder_intersector ci;

    // Unbound local frame test
    unbound ub;
    auto hit_unbound = ci.intersect(shifted_cylinder, point3{3., 2., 5.}, vector3{1., 0., 0.}, ctx, ub, cylinder);
    ASSERT_TRUE(hit_unbound._status == intersection_status::e_inside);
    ASSERT_NEAR(hit_unbound._point3[0], 7., epsilon);
    ASSERT_NEAR(hit_unbound._point3[1], 2., epsilon);
    ASSERT_NEAR(hit_unbound._point3[2], 5., epsilon);
    ASSERT_TRUE(hit_unbound._point2 == std::nullopt);

    // The same but bound
    plugin::cylindrical2 cylindrical2;
    auto hit_bound = ci.intersect(shifted_cylinder, point3{3., 2., 5.}, vector3{1., 0., 0.}, ctx, cylindrical2, cylinder);
    ASSERT_TRUE(hit_bound._status == intersection_status::e_inside);
    ASSERT_NEAR(hit_bound._point3[0], 7., epsilon);
    ASSERT_NEAR(hit_bound._point3[1], 2., epsilon);
    ASSERT_NEAR(hit_bound._point3[2], 5., epsilon);
    ASSERT_TRUE(hit_bound._point2 != std::nullopt);
    ASSERT_NEAR(hit_bound._point2.value()[0], 0., isclose);
    ASSERT_NEAR(hit_bound._point2.value()[1], -5., isclose);
}

// This defines the local frame test suite
TEST(plugin, concentric_cylinders)
{
    context ctx;
    using cylinder_surface = surface<transform3, int, int>;

    // Create a shifted plane
    scalar r = 4.;
    scalar hz = 10.;
    transform3 identity(vector3{0., 0., 0.}, ctx);
    cylinder_surface plain(std::move(identity), 1, 1);
    cylinder3<scalar> cylinder = {r, hz};
    single3<scalar, 2> halflength = {hz};
    cylinder_intersector ci;
    concentric_cylinder_intersector cci;

    point3 ori = {1., 0.5, 1.};
    point3 dir = vector::normalize(vector3{1., 1., 1.});

    // The same but bound
    plugin::cylindrical2 cylindrical2;
    auto hit_cylinrical = ci.intersect(plain, ori, dir, ctx, cylindrical2, cylinder);
    auto hit_cocylinrical = cci.intersect(plain, r, ori, dir, ctx, cylindrical2, halflength);

    ASSERT_TRUE(hit_cylinrical._status == intersection_status::e_inside);
    ASSERT_TRUE(hit_cocylinrical._status == intersection_status::e_inside);

    ASSERT_NEAR(getter::perp(hit_cylinrical._point3), r, isclose);
    ASSERT_NEAR(getter::perp(hit_cocylinrical._point3), r, isclose);

    ASSERT_NEAR(hit_cylinrical._point3[0], hit_cocylinrical._point3[0], isclose);
    ASSERT_NEAR(hit_cylinrical._point3[1], hit_cocylinrical._point3[1], isclose);
    ASSERT_NEAR(hit_cylinrical._point3[2], hit_cocylinrical._point3[2], isclose);
    ASSERT_TRUE(hit_cylinrical._point2 != std::nullopt);
    ASSERT_TRUE(hit_cocylinrical._point2 != std::nullopt);
    ASSERT_NEAR(hit_cylinrical._point2.value()[0], hit_cocylinrical._point2.value()[0], isclose);
    ASSERT_NEAR(hit_cylinrical._point2.value()[1], hit_cocylinrical._point2.value()[1], isclose);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
