/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "tools/navigator.hpp"
#include "tests/common/test_detector.hpp"

#include <fstream>
#include <cmath>
#include <climits>

#include <gtest/gtest.h>

unsigned int theta_steps = 10;
unsigned int phi_steps = 10000;
bool stream_file = true;

auto d = createDetector();

TEST(__plugin, intersect_all_cylindrical_detector)
{
    navigator<cdetector> n;

    cdetector::surface_intersection sfi;
    const auto &surfaces = d.surfaces();
    const auto &surface_transforms = d.surface_transforms();
    const auto &stm = d.surface_types();
    const auto &msks = d.surface_masks();
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
                sfi.index = s;
                const auto &surface = surfaces[s];
                const auto &transform = surface_transforms[surface.transform()];

                if (n.update_intersection(sfi, links, track, transform, surface, stm, msks) and stream_file)
                {
                    hit_out << sfi.point3[0] << "," << sfi.point3[1] << "," << sfi.point3[2] << "\n";
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

    navigator<cdetector> n;
    navigator<cdetector>::navigation_state navigation;
    navigation.detector = d;
    navigation.volume_index = 0; // this will have to be found out by the detector... octotree?

    // Starting form 0,0,0 with direction (1,1,1).normalized
    track<transform3> track;
    track.pos = {0., 0., 0.};
    track.dir = {1. / std::sqrt(3), 1. / std::sqrt(3), 1. / std::sqrt(3)};
    // Target to next
    ASSERT_TRUE(n.target(navigation, track) == navigator<cdetector>::navigation_status::e_towards_surface);
    track.pos = track.pos + navigation.distance_to_next * track.dir;
    ASSERT_TRUE(n.status(navigation, track, true) == navigator<cdetector>::navigation_status::e_on_surface);
    // Let's target again
    ASSERT_TRUE(n.target(navigation, track) == navigator<cdetector>::navigation_status::e_towards_portal);
    track.pos = track.pos + navigation.distance_to_next * track.dir;
    ASSERT_TRUE(n.status(navigation, track, true) == navigator<cdetector>::navigation_status::e_on_portal);
    // Target to next
    ASSERT_TRUE(n.target(navigation, track) == navigator<cdetector>::navigation_status::e_towards_surface);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
