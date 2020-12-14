/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/surface.hpp"
#include "geometry/cylindrical_detector.hpp"
#include "masks/cylinder3.hpp"
#include "masks/ring2.hpp"

#include <cmath>
#include <climits>

#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

// Three-dimensional definitions
using transform3 = __plugin::transform3;
using vector3 = __plugin::transform3::vector3;
using point3 = __plugin::transform3::point3;
using context = __plugin::transform3::context;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();
constexpr scalar isclose = 1e-5;

// This defines the local frame test suite
TEST(__plugin, cylindrical_detector)
{

    using portal_range = darray<int, 2>;
    using portal_links = dtuple<unsigned int, portal_range>;
    using surface = surface<transform3, portal_links, unsigned int>;
    using detector = cylindrical_detector<transform3>;

    context ctx;
    detector d;

    scalar bp_radius = 32.;
    scalar bp_length = 1000.;
    scalar p_radius = 120.;
    scalar p_barrel = 400.;

    detector::portal_cylinder_mask bp_cylinder_ecn = {{bp_radius, -0.5 * bp_length, -0.5 * p_barrel}, {0, 1}};
    detector::portal_cylinder_mask bp_cylinder_b = {{bp_radius, -0.5 * p_barrel, 0.5 * p_barrel}, {0, 2}};
    detector::portal_cylinder_mask bp_cylinder_ecp = {{bp_radius, 0.5 * p_barrel, 0.5 * bp_length}, {0, 3}};
    detector::portal_disc_mask bp_n_disc = {{0., 23.}, {-1, 0}};
    detector::portal_disc_mask bp_p_disc = {{0., 23.}, {0, -1}};

    unsigned int sfcoutner = 0;

    detector::volume bp;
    bp.volume_index = 0;
    dvector<detector::portal_cylinder_mask> bp_p_portals = {bp_cylinder_ecn, bp_cylinder_b, bp_cylinder_ecp};
    auto bp_b_portal_links = d.add_portals<detector::portal_cylinder_mask>(bp_p_portals);
    detector::portal_surface bp_b(std::move(transform3{}), std::move(bp_b_portal_links), std::move(false));
    detector::guaranteed_index start_index = d.add_portal_surface(std::move(bp_b));
    auto bp_ecn_portal_links = d.add_portals<detector::portal_disc_mask>({bp_n_disc});
    detector::portal_surface bp_ecn(std::move(transform3(vector3(0., 0., 0.5 * bp_length), ctx)), std::move(bp_ecn_portal_links), std::move(false));
    detector::guaranteed_index current_index = d.add_portal_surface(std::move(bp_ecn));
    auto bp_ecp_portal_links = d.add_portals<detector::portal_disc_mask>({bp_p_disc});
    detector::portal_surface bp_ecp(std::move(transform3(vector3(0., 0., 0.5 * bp_length), ctx)), std::move(bp_ecp_portal_links), std::move(false));
    current_index = d.add_portal_surface(std::move(bp_ecp));
    bp.portal_surface_range = {start_index, current_index};
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
