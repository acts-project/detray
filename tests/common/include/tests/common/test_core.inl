/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/surface.hpp"
#include "core/intersection.hpp"

#include <cmath>
#include <climits>

#include <gtest/gtest.h>

/// @note plugin has to be defined with a preprocessor command

using namespace detray;

using point2 = plugin::cartesian2::point2;

// Three-dimensional definitions
using transform3 = plugin::transform3;
using vector3 = plugin::transform3::vector3;
using point3 = plugin::transform3::point3;
using context = plugin::transform3::context;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();

// This tests the construction of a surface
TEST(plugin, surface)
{
    context ctx;

    // Preparatioon work, create a transform
    vector3 z = vector::normalize(vector3(3., 2., 1.));
    vector3 x = vector::normalize(vector3(2., -3., 0.));
    vector3 y = vector::cross(z, x);
    point3 t(2., 3., 4.);
    transform3 trf(t, z, x, ctx);

    surface s(std::move(trf),-1,-1);
}

// This tests the construction of a intresection
TEST(plugin, intersection)
{
    using intersection = intersection<scalar, point3, point2>;

    intersection i0 = {2., point3(0.3, 0.5, 0.7), point2(0.2, 0.4), intersection_status::e_hit};

    intersection i1 = {1.7, point3(0.2, 0.3, 0.), point2(0.2, 0.4), intersection_status::e_inside};

    intersection invalid;
    ASSERT_TRUE(invalid._status == intersection_status::e_missed);

    std::vector<intersection> intersections = {invalid, i0, i1};
    std::sort(intersections.begin(), intersections.end());

    ASSERT_NEAR(intersections[0]._path, 1.7, epsilon);
    ASSERT_NEAR(intersections[1]._path, 2, epsilon);
    ASSERT_TRUE(std::isinf(intersections[2]._path));
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
