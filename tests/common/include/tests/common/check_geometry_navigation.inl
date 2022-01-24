/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/track.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/tools/line_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/ray_gun.hpp"

using namespace detray;

/** A navigation inspector that relays information about the encountered
 *  objects the way we need them to compare with the ray
 */
template <int navigation_status = 0,
          template <typename...> class vector_t = dvector>
struct object_tracer {

    // record all object id the navigator encounters
    vector_t<intersection> object_trace = {};

    template <typename state_type>
    auto operator()(state_type &state, std::string & /*message*/) {
        // Record the candidate of an encountered object
        if (state.status() == navigation_status) {
            object_trace.push_back(std::move(*(state.current())));
        }
    }

    auto operator[](std::size_t i) { return object_trace[i]; }
};

/** A navigation inspector that prints information about the current navigation
 * state. Meant for debugging.
 */
struct print_inspector {

    // Debug output if an error in the trace is discovered
    std::stringstream debug_stream;

    template <typename state_type>
    auto operator()(state_type &state, std::string &message) {
        debug_stream << message << std::endl;

        debug_stream << "Volume\t\t\t\t\t\t" << state.volume() << std::endl;
        debug_stream << "surface kernel size\t\t" << state.candidates().size()
                     << std::endl;

        debug_stream << "Surface candidates: " << std::endl;
        for (const auto &sf_cand : state.candidates()) {
            debug_stream << sf_cand.to_string();
        }
        if (not state.candidates().empty()) {
            debug_stream << "=> next: ";
            if (state.is_exhausted()) {
                debug_stream << "exhausted" << std::endl;
            } else {
                debug_stream << " -> " << state.next()->index << std::endl;
            }
        }

        switch (static_cast<int>(state.status())) {
            case -3:
                debug_stream << "status\t\t\t\t\ton_target" << std::endl;
                break;
            case -2:
                debug_stream << "status\t\t\t\t\tabort" << std::endl;
                break;
            case -1:
                debug_stream << "status\t\t\t\t\tunknowm" << std::endl;
                break;
            case 0:
                debug_stream << "status\t\t\t\t\ttowards_surface" << std::endl;
                break;
            case 1:
                debug_stream << "status\t\t\t\t\ton_surface" << std::endl;
                break;
            case 2:
                debug_stream << "status\t\t\t\t\ttowards_portal" << std::endl;
                break;
            case 3:
                debug_stream << "status\t\t\t\t\ton_portal" << std::endl;
                break;
        };
        debug_stream << "current object\t\t" << state.on_object() << std::endl;
        debug_stream << "distance to next\t";
        if (std::abs(state()) < state.tolerance()) {
            debug_stream << "on obj (within tol)" << std::endl;
        } else {
            debug_stream << state() << std::endl;
        }
        switch (state.trust_level()) {
            case 0:
                debug_stream << "trust\t\t\t\t\tno_trust" << std::endl;
                break;
            case 1:
                debug_stream << "trust\t\t\t\t\tfair_trust" << std::endl;
                break;
            case 3:
                debug_stream << "trust\t\t\t\t\thigh_trust" << std::endl;
                break;
            case 4:
                debug_stream << "trust\t\t\t\t\tfull_trust" << std::endl;
                break;
        };
        debug_stream << std::endl;
    }

    std::string to_string() { return debug_stream.str(); }
};

/** A navigation inspector that aggregates a number of different inspectors.*/
template <typename... Inspectors>
struct aggregate_inspector {

    using inspector_tuple_t = std::tuple<Inspectors...>;
    inspector_tuple_t _inspectors{};

    template <unsigned int current_id = 0, typename state_type>
    auto operator()(state_type &state, std::string &message) {
        // Call inspector
        std::get<current_id>(_inspectors)(state, message);

        // Next mask type
        if constexpr (current_id <
                      std::tuple_size<inspector_tuple_t>::value - 1) {
            return operator()<current_id + 1>(state, message);
        }
    }

    template <typename inspector_t>
    decltype(auto) get() {
        return std::get<inspector_t>(_inspectors);
    }
};

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, geometry_discovery) {
    // vecmem::host_memory_resource host_mr;
    // auto [d, name_map] = read_from_csv(tml_files, host_mr);

    vecmem::host_memory_resource host_mr;
    auto d = create_toy_geometry(host_mr);

    // Create the navigator
    using detray_context = decltype(d)::context;
    using detray_track = track<detray_context>;
    using detray_inspector =
        aggregate_inspector<object_tracer<1>, print_inspector>;
    using detray_navigator = navigator<decltype(d), detray_inspector>;
    using detray_stepper = line_stepper<detray_track>;

    detray_navigator n(d);
    detray_stepper s;

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

            // Now follow that ray and check, if we find the same
            // volumes and distances along the way
            track<detray_context> ray;
            ray.pos = ori;
            ray.dir = dir;
            ray.ctx = detray_context{};
            ray.momentum = 100.;
            ray.overstep_tolerance = -1e-4;

            std::vector<std::pair<dindex, intersection>> intersection_trace;

            shoot_ray(d, ray, intersection_trace);

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

            auto &obj_tracer =
                n_state.inspector().template get<object_tracer<1>>();
            auto &debug_printer =
                n_state.inspector().template get<print_inspector>();

            std::stringstream debug_stream;
            for (std::size_t intr_idx = 0; intr_idx < intersection_trace.size();
                 ++intr_idx) {
                debug_stream
                    << "-------Intersection trace\n"
                    << "ray gun: "
                    << "\tsf id: " << intersection_trace[intr_idx].first << ", "
                    << intersection_trace[intr_idx].second.to_string();
                debug_stream << "navig.: " << obj_tracer[intr_idx].to_string();
            }

            // Compare intersection records
            EXPECT_EQ(obj_tracer.object_trace.size(),
                      intersection_trace.size());
            // Check every single recorded intersection
            for (std::size_t intr_idx = 0; intr_idx < intersection_trace.size();
                 ++intr_idx) {
                if (obj_tracer[intr_idx].index !=
                    intersection_trace[intr_idx].first) {
                    // Intersection record at portal bound might be flipped
                    if (obj_tracer[intr_idx].index ==
                            intersection_trace[intr_idx + 1].first and
                        obj_tracer[intr_idx + 1].index ==
                            intersection_trace[intr_idx].first) {
                        // Have already checked the next record
                        ++intr_idx;
                        continue;
                    }
                }
                EXPECT_EQ(obj_tracer[intr_idx].index,
                          intersection_trace[intr_idx].first)
                    << debug_printer.to_string() << debug_stream.str();
            }
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
