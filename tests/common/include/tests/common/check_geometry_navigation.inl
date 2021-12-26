/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/track.hpp"
#include "detray/tools/line_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "tests/common/tools/ray_gun.hpp"
//#include "tests/common/tools/read_geometry.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

using namespace detray;

/** A navigation inspector that relays information about the encountered
 *  objects the way we need them to compare with the ray
 */
template <int navigation_status = 0,
          template <typename...> class vector_t = dvector>
struct object_tracer {

    // record all object id the navigator encounters
    vector_t<dindex> object_trace = {};

    template <typename state_type>
    auto operator()(state_type &state, std::string & /*message*/) {
        // Record and objects id, when you are certain to have encountered it
        if (state.status() == navigation_status) {
            object_trace.push_back(state.on_object());
        }
    }
};

/** A navigation inspector that prints information about the current navigation
 * state. Meant for debugging.
 */
struct print_inspector {

    template <typename state_type>
    auto operator()(state_type &state, std::string &message) {

        std::cout << message << std::endl;

        std::cout << "Volume\t\t\t\t\t" << state.volume() << std::endl;
        std::cout << "surface kernel size\t\t" << state.kernel().size()
                  << std::endl;

        std::cout << "Surface candidates: " << std::endl;
        for (const auto &sf_cand : state.candidates()) {
            std::cout << "-> " << sf_cand.path << " (" << sf_cand.index
                      << ", links to " << sf_cand.link << ")" << std::endl;
        }
        if (not state.kernel().empty()) {
            std::cout << "=> next: ";
            if (state.is_exhausted()) {
                std::cout << "exhausted" << std::endl;
            } else {
                std::cout << state.next()->index << std::endl;
            }
        }

        switch (state.status()) {
            case -3:
                std::cout << "status\t\t\t\ton_target" << std::endl;
                break;
            case -2:
                std::cout << "status\t\t\t\tabort" << std::endl;
                break;
            case -1:
                std::cout << "status\t\t\t\tunknowm" << std::endl;
                break;
            case 0:
                std::cout << "status\t\t\t\ttowards_surface" << std::endl;
                break;
            case 1:
                std::cout << "status\t\t\t\ton_surface" << std::endl;
                break;
            case 2:
                std::cout << "status\t\t\t\ttowards_portal" << std::endl;
                break;
            case 3:
                std::cout << "status\t\t\t\ton_portal" << std::endl;
                break;
        };
        std::cout << "current object\t\t" << state.on_object() << std::endl;
        std::cout << "distance to next\t";
        if (std::abs(state()) < state.tolerance()) {
            std::cout << "on obj (within tol)" << std::endl;
        } else {
            std::cout << state() << std::endl;
        }
        switch (state.trust_level()) {
            case 0:
                std::cout << "trust\t\t\t\tno_trust" << std::endl;
                break;
            case 1:
                std::cout << "trust\t\t\t\tfair_trust" << std::endl;
                break;
            case 3:
                std::cout << "trust\t\t\t\thigh_trust" << std::endl;
                break;
            case 4:
                std::cout << "trust\t\t\t\tfull_trust" << std::endl;
                break;
        };
        std::cout << std::endl;
    }
};

// vecmem::host_memory_resource host_mr;
// auto [d, name_map] = read_from_csv(tml_files, host_mr);

vecmem::host_memory_resource host_mr;
auto d = create_toy_geometry(host_mr);

// Create the navigator
using detray_context = decltype(d)::context;
using detray_track = track<detray_context>;
using detray_print_inspector = print_inspector;
using detray_inspector = object_tracer<1>;
using detray_navigator = navigator<decltype(d), detray_inspector>;
using detray_stepper = line_stepper<detray_track>;

detray_navigator n(d);
detray_stepper s;

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, geometry_discovery) {

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

            const auto intersection_trace = shoot_ray(d, ori, dir);

            // Now follow that ray and check, if we find the same
            // volumes and distances along the way
            track<detray_context> ray;
            ray.pos = ori;
            ray.dir = dir;
            ray.ctx = detray_context{};
            ray.momentum = 100.;
            ray.overstep_tolerance = -1e-4;

            detray_stepper::state s_state(ray);
            detray_navigator::state n_state;
            // Always start a new ray at detector origin
            n_state.set_volume(0u);

            bool heartbeat = n.status(n_state, ray);
            // Run while there is a heartbeat
            while (heartbeat) {
                // (Re-)target
                heartbeat &= n.target(n_state, s_state());
                // Take the step
                heartbeat &= s.step(s_state, n_state());
                // And check the status
                heartbeat &= n.status(n_state, s_state());
            }
            // Compare intersection records
            if constexpr (std::is_same_v<detray_inspector, object_tracer<1>>) {
                EXPECT_EQ(n_state.inspector().object_trace.size(),
                          intersection_trace.size());
                // Check every single recorded intersection
                for (std::size_t intr_idx = 0;
                     intr_idx < intersection_trace.size(); ++intr_idx) {
                    if (n_state.inspector().object_trace[intr_idx] !=
                              intersection_trace[intr_idx].first) {
                        // Intersection record at portal bound might be flipped
                        if (n_state.inspector().object_trace[intr_idx] == intersection_trace[intr_idx + 1].first and n_state.inspector().object_trace[intr_idx + 1] == intersection_trace[intr_idx].first) {
                            // Have already checked the next record
                            ++intr_idx;
                            continue;
                        }
                    }
                    EXPECT_EQ(n_state.inspector().object_trace[intr_idx],
                              intersection_trace[intr_idx].first);
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
