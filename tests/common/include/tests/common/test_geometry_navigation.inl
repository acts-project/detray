/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "core/detector.hpp"
#include "core/track.hpp"
#include "io/csv_io.hpp"
#include "tools/line_stepper.hpp"
#include "tools/navigator.hpp"
#include "utils/ray_gun.hpp"

using namespace detray;

/** Read the detector from file */
auto read_detector() {
    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr) {
        throw std::ios_base::failure(
            "Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string name = "tml";
    std::string surface_file = data_directory + "tml.csv";
    std::string layer_volume_file = data_directory + "tml-layer-volumes.csv";
    std::string surface_grid_file = data_directory + "tml-surface-grids.csv";
    std::string surface_grid_entries_file = "";

    /*std::string name = "odd";
    std::string surface_file = data_directory + std::string("odd.csv");
    std::string layer_volume_file =
        data_directory + std::string("odd-layer-volumes.csv");
    std::string surface_grid_file =
        data_directory + std::string("odd-surface-grids.csv");
    std::string surface_grid_entries_file = "";*/

    return detector_from_csv<>(name, surface_file, layer_volume_file,
                               surface_grid_file, surface_grid_entries_file);
};

/** A navigation inspector that relays information about the encountered
 *  portals the way we need them to compare with the ray
 */
struct portal_inspector {
    template <typename state_type>
    auto operator()(const state_type &state) {
        // on portal
        /*if (state.status == 3) {
            std::cout << "Hit portal" << std::endl;
        }*/
        std::cout << state.status << std::endl;
        std::cout << state.volume_index << std::endl;
        std::cout << state.distance_to_next << std::endl;
    }
};

auto d = read_detector();

// Create the navigator
using detray_context = decltype(d)::context;
using detray_track = track<detray_context>;
using detray_navigator = navigator<decltype(d), portal_inspector>;
using detray_stepper = line_stepper<detray_track>;

detray_navigator n(std::move(d));
detray_stepper s;

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, ray_scan) {

    /*unsigned int theta_steps = 100;
    unsigned int phi_steps = 100;
    const unsigned int itest = 10000;*/
    unsigned int theta_steps = 1;
    unsigned int phi_steps = 1;

    const point3 ori{0., 0., 0.};
    dindex start_index = d.volume_by_pos(ori).index();

    // Loops of theta values ]0,pi[
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        // Loops of phi values [-pi, pi]
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {
            // The direction
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            scalar sin_phi = std::sin(phi);
            scalar cos_phi = std::cos(phi);
            const point3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                             cos_theta};

            const auto volume_record = shoot_ray(d, ori, dir);

            // Now follow that ray and check, if we find the same
            // volumes and distances along the way
            track<detray_context> traj;
            traj.pos = ori;
            traj.dir = dir;
            traj.ctx = detray_context{};
            traj.momentum = 100.;
            traj.overstep_tolerance = 0.;

            detray_stepper::state s_state(traj);
            detray_navigator::state n_state;

            bool heartbeat = n.status(n_state, s_state());
            // Run while there is a heartbeat
            while (heartbeat) {
                // (Re-)target
                heartbeat &= n.target(n_state, s_state());
                // Take the step
                heartbeat &= s.step(s_state, n_state());
                // And check the status
                heartbeat &= n.status(n_state, s_state());
            }
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
