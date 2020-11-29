
#include "core/surface.hpp"
#include "core/intersection.hpp"
#include "masks/rectangle2.hpp"
#include "tools/planar_intersector.hpp"

#include <cmath>
#include <climits>

/// @note plugin has to be defined with a preprocessor command
using namespace detray;

// Two-dimensional definitions
using point2cart = plugin::cartesian2::point2;
plugin::cartesian2 cartesian2;

// Three-dimensional definitions
using transform3 = plugin::transform3;
using vector3 = plugin::transform3::vector3;
using point3 = plugin::transform3::point3;
using context = plugin::transform3::context;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();
constexpr scalar isclose = 1e-5;


// This defines the local frame test suite
TEST(plugin, translated_plane)
{
    context ctx;
    using plane_surface = surface<transform3,int>;

    // Create a shifted plane
    transform3 shifted(vector3{3., 2., 10.}, ctx);
    plane_surface shifted_plane(std::move(shifted), 1);

    planar_intersector pi;

    auto hit0 = pi.intersect(shifted_plane, point3{2.,1.,0.}, vector3{0.,0.,1.}, ctx);
    ASSERT_NEAR(hit0._point3[0], 2., epsilon);
    ASSERT_NEAR(hit0._point3[1], 1., epsilon);
    ASSERT_NEAR(hit0._point3[2], 10., epsilon);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
