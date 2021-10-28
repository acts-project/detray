/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vecmem/memory/host_memory_resource.hpp>

#include "core/detector.hpp"
#include "core/track.hpp"
#include "io/csv_io.hpp"
#include "tests/common/read_geometry.hpp"
#include "tools/line_stepper.hpp"
#include "tools/navigator.hpp"
#include "utils/ray_gun.hpp"

using namespace detray;

vecmem::host_memory_resource host_mr;

/** A navigation inspector that relays information about the encountered
 *  portals the way we need them to compare with the ray
 */
struct portal_inspector {
    template <typename state_type>
    auto operator()(const state_type &state) {}
};

/** A navigation inspector that prints information about the current navigation
 * state. Meant for debugging.
 */
struct print_inspector {
    template <typename state_type>
    auto operator()(const state_type &state) {

        std::cout << "Volume\t\t\t" << state.volume_index << std::endl;
        std::cout << "surface kernel size\t\t\t" << state.surface_kernel.size()
                  << std::endl;
        std::cout << "portal kernel size\t\t\t" << state.portal_kernel.size()
                  << std::endl;

        std::cout << "Surface candidates: " << std::endl;
        for (const auto &sf_cand : state.surface_kernel.candidates) {
            std::cout << "-> " << sf_cand.path << std::endl;
        }
        if (not state.surface_kernel.empty())
            std::cout << "=> next: " << state.surface_kernel.next->index
                      << std::endl;

        std::cout << "Portal candidates: " << std::endl;
        for (const auto &pt_cand : state.portal_kernel.candidates) {
            std::cout << "-> " << pt_cand.path << std::endl;
        }
        if (not state.portal_kernel.empty())
            std::cout << "=> next: " << state.portal_kernel.next->index
                      << std::endl;

        switch (state.status) {
            case -3:
                std::cout << "status\t\t\ton_target" << std::endl;
                break;
            case -2:
                std::cout << "status\t\t\tabort" << std::endl;
                break;
            case -1:
                std::cout << "status\t\t\tunknowm" << std::endl;
                break;
            case 0:
                std::cout << "status\t\t\ttowards_surface" << std::endl;
                break;
            case 1:
                std::cout << "status\t\t\ton_surface" << std::endl;
                break;
            case 2:
                std::cout << "status\t\t\ttowards_portal" << std::endl;
                break;
            case 3:
                std::cout << "status\t\t\ton_portal" << std::endl;
                break;
        };
        std::cout << "current object\t\t" << state.current_index << std::endl;
        std::cout << "distance to next\t" << state.distance_to_next
                  << std::endl;
        switch (state.trust_level) {
            case 0:
                std::cout << "trust\t\t\tno_trust" << std::endl;
                break;
            case 1:
                std::cout << "trust\t\t\tfair_trust" << std::endl;
                break;
            case 3:
                std::cout << "trust\t\t\thigh_trust" << std::endl;
                break;
            case 4:
                std::cout << "trust\t\t\tfull_trust" << std::endl;
                break;
        };
        std::cout << std::endl;
    }
};

auto [d, name_map] = read_from_csv(tml_files, host_mr);

// Create the navigator
using detray_context = decltype(d)::context;
using detray_track = track<detray_context>;
using detray_navigator = navigator<decltype(d), print_inspector>;
using detray_stepper = line_stepper<detray_track>;

detray_navigator n(std::move(d));
detray_stepper s;

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, ray_scan) {

    unsigned int theta_steps = 100;
    unsigned int phi_steps = 100;

    const point3 ori{0., 0., 0.};
    // dindex start_index = n.detector.volume_by_pos(ori).index();

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

            const auto volume_record = shoot_ray(n.detector, ori, dir);

            // Now follow that ray and check, if we find the same
            // volumes and distances along the way
            track<detray_context> traj;
            traj.pos = ori;
            traj.dir = dir;
            traj.ctx = detray_context{};
            traj.momentum = 100.;
            traj.overstep_tolerance = -1e-4;

            detray_stepper::state s_state(traj);
            detray_navigator::state n_state;

            bool heartbeat = n.status(n_state, s_state());
            // Run while there is a heartbeat
            // while (heartbeat) {
            for (size_t n_steps = 0; n_steps < 10; n_steps++) {
                // (Re-)target
                heartbeat &= n.target(n_state, s_state());
                // Take the step
                heartbeat &= s.step(s_state, n_state());
                // And check the status
                heartbeat &= n.status(n_state, s_state());
            }
            // EXPECT_EQ(n_state.current_index, volume_record.back().first);
            // EXPECT_EQ(getter::norm(traj.pos),
            // volume_record.back().second.path);
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
