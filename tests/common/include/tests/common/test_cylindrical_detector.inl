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
#include "tests/common/test_surfaces.hpp"

#include <fstream>
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

using detector = cylindrical_detector<transform3>;
detector d;
context ctx;

// This test constructs a dimple cylindrical detector
TEST(__plugin, cylindrical_detector)
{
    // An inner volume: call it pb
    scalar bp_radius = 29.;
    scalar bp_length = 1000.;
    scalar px_barrel = 600.;
    scalar px_endcap = 0.5 * (bp_length - px_barrel);

    detector::volume &bp = d.new_volume("bp", {0., bp_radius, -0.5 * bp_length, 0.5 * bp_length});
    detector::portal_cylinder_mask bp_c_ecn = {{bp_radius, -0.5 * bp_length, -0.5 * px_barrel}, {0, 1}};
    detector::portal_cylinder_mask bp_c_b = {{bp_radius, -0.5 * px_barrel, 0.5 * px_barrel}, {0, 2}};
    detector::portal_cylinder_mask bp_c_ecp = {{bp_radius, 0.5 * px_barrel, 0.5 * bp_length}, {0, 3}};
    detector::portal_disc_mask bp_n_disc = {{0., bp_radius}, {-1, 0}};
    detector::portal_disc_mask bp_p_disc = {{0., bp_radius}, {0, -1}};
    dvector<detector::portal_cylinder_mask> bp_c_portals = {bp_c_ecn, bp_c_b, bp_c_ecp};
    d.add_portal_surface<detector::portal_cylinder_mask>(std::move(transform3()), bp_c_portals, bp);
    d.add_portal_surface<detector::portal_disc_mask>(std::move(transform3(vector3(0., 0., -0.5 * bp_length), ctx)), {bp_n_disc}, bp);
    d.add_portal_surface<detector::portal_disc_mask>(std::move(transform3(vector3(0., 0., 0.5 * bp_length), ctx)), {bp_p_disc}, bp);
    // Insert an actual beam pipe
    darray<scalar, 3> bpm_values = {25., -0.5 * bp_length + 1., 0.5 * bp_length - 1.};
    d.add_surfaces<detector::cylinder_mask>({transform3()}, bpm_values, bp);

    // A wrapping volume, let's call it px
    scalar px_inner_radius = bp_radius;
    scalar px_outer_radius = 55.;

    detector::volume &px_ecn = d.new_volume("px_ecn", {px_inner_radius, px_outer_radius, -bp_length, -px_barrel});
    detector::portal_cylinder_mask px_ecn_inner = {{px_inner_radius, -bp_length, -px_barrel}, {0, 1}};
    detector::portal_cylinder_mask px_ecn_outer = {{px_outer_radius, -bp_length, -px_barrel}, {1, -1}};
    detector::portal_disc_mask px_ecn_ecn = {{px_inner_radius, px_outer_radius}, {-1, 1}};
    detector::portal_disc_mask px_ecn_ecp = {{px_inner_radius, px_outer_radius}, {1, 2}};
    d.add_portal_surface<detector::portal_cylinder_mask>(std::move(transform3()), {px_ecn_inner, px_ecn_outer}, px_ecn);
    d.add_portal_surface<detector::portal_disc_mask>(std::move(transform3(vector3(0., 0., -0.5 * bp_length), ctx)), {px_ecn_ecn}, px_ecn);
    d.add_portal_surface<detector::portal_disc_mask>(std::move(transform3(vector3(0., 0., -0.5 * px_barrel), ctx)), {px_ecn_ecp}, px_ecn);
    auto endcap_n = endcap_description(29, 50, -px_barrel - 25., 1., 6, 0.2);
    d.add_surfaces<detector::trapezoid_mask>(std::get<dvector<transform3>>(endcap_n),
                                             std::get<detector::trapezoid_mask::mask_values>(endcap_n),
                                             px_ecn);

    detector::volume &px_b = d.new_volume("px_b", {px_inner_radius, px_outer_radius, -px_barrel, px_barrel});
    detector::portal_cylinder_mask px_b_inner = {{px_inner_radius, -px_barrel, px_barrel}, {0, 2}};
    detector::portal_cylinder_mask px_b_outer = {{px_outer_radius, -px_barrel, px_barrel}, {2, -1}};
    detector::portal_disc_mask px_b_ecn = {{px_inner_radius, px_outer_radius}, {1, 2}};
    detector::portal_disc_mask px_b_ecp = {{px_inner_radius, px_outer_radius}, {2, 3}};
    d.add_portal_surface<detector::portal_cylinder_mask>(std::move(transform3()), {px_b_inner, px_b_outer}, px_b);
    d.add_portal_surface<detector::portal_disc_mask>(std::move(transform3(vector3(0., 0., -0.5 * px_barrel), ctx)), {px_b_ecn}, px_b);
    d.add_portal_surface<detector::portal_disc_mask>(std::move(transform3(vector3(0., 0., +0.5 * px_barrel), ctx)), {px_b_ecp}, px_b);
    auto barrel = barrel_description(33, 0.25, 12, 0.12, 0.25, 600, 2., 7);
    d.add_surfaces<detector::rectangle_mask>(std::get<dvector<transform3>>(barrel),
                                             std::get<detector::rectangle_mask::mask_values>(barrel),
                                             px_b);

    detector::volume &px_ecp = d.new_volume("px_ecp", {px_inner_radius, px_outer_radius, px_barrel, bp_length});
    detector::portal_cylinder_mask px_ecp_inner = {{px_inner_radius, px_barrel, bp_length}, {0, 3}};
    detector::portal_cylinder_mask px_ecp_outer = {{px_outer_radius, px_barrel, bp_length}, {3, -1}};
    detector::portal_disc_mask px_ecp_ecn = {{px_inner_radius, px_outer_radius}, {2, 3}};
    detector::portal_disc_mask px_ecp_ecp = {{px_inner_radius, px_outer_radius}, {3, -1}};
    d.add_portal_surface<detector::portal_cylinder_mask>(std::move(transform3()), {px_ecp_inner, px_ecp_outer}, px_ecp);
    d.add_portal_surface<detector::portal_disc_mask>(std::move(transform3(vector3(0., 0., 0.5 * px_barrel), ctx)), {px_ecp_ecn}, px_ecp);
    d.add_portal_surface<detector::portal_disc_mask>(std::move(transform3(vector3(0., 0., 0.5 * bp_length), ctx)), {px_ecp_ecp}, px_ecp);
    auto endcap_p = endcap_description(29, 50, px_barrel + 25., 1., 6, 0.2);
    d.add_surfaces<detector::trapezoid_mask>(std::get<dvector<transform3>>(endcap_p),
                                             std::get<detector::trapezoid_mask::mask_values>(endcap_p),
                                             px_ecp);
}

unsigned int theta_steps = 10;
unsigned int phi_steps = 10000;
bool stream_file = true;

TEST(__plugin, intersect_random_cylindrical_detector)
{
    detector::surface_intersection sfi;
    auto &surfaces = d._surfaces;
    const auto &stm = d._surface_types;
    const auto &msks = d._surface_masks;
    bool links = false;

    std::ofstream hit_out;
    if (stream_file)
    {
        hit_out.open("three_layers.csv");
    }

    track<transform3> track;
    track.pos = {0., 0., 0.};

    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta)
    {
        scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);

        // Loops of phi values
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi)
        {
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            double sin_phi = std::sin(phi);
            double cos_phi = std::cos(phi);

            track.dir = {cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};

            for (unsigned int s = 0; s < surfaces.size(); ++s)
            {
                sfi._index = s;
                if (d.update_intersection(sfi, links, track, surfaces, stm, msks) and stream_file)
                {
                    hit_out << sfi._point3[0] << "," << sfi._point3[1] << "," << sfi._point3[2] << "\n";
                }
            }
        }
    }
    if (stream_file)
    {
        hit_out.close();
    }
}

// This test navigates through a cylindrical detector
TEST(__plugin, navigate_cylindrical_detector)
{
    detector::navigation_state navigation;
    navigation.volume_index = 0; // this will have to be found out by the detector... octotree?

    // Starting form 0,0,0 with direction (1,1,1).normalized
    track<transform3> track;
    track.pos = {0., 0., 0.};
    track.dir = {1. / std::sqrt(3), 1. / std::sqrt(3), 1. / std::sqrt(3)};
    // Target to next
    ASSERT_TRUE(d.target(navigation, track) == detector::navigation_status::e_towards_surface);
    track.pos = track.pos + navigation.distance_to_next * track.dir;
    ASSERT_TRUE(d.status(navigation, track, true) == detector::navigation_status::e_on_surface);
    // Let's target again
    ASSERT_TRUE(d.target(navigation, track) == detector::navigation_status::e_towards_portal);
    track.pos = track.pos + navigation.distance_to_next * track.dir;
    ASSERT_TRUE(d.status(navigation, track, true) == detector::navigation_status::e_on_portal);
    // Target to next
    ASSERT_TRUE(d.target(navigation, track) == detector::navigation_status::e_towards_surface);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
