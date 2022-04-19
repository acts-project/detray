/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <sstream>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/track.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/inspectors.hpp"
#include "tests/common/tools/ray_gun.hpp"

using namespace detray;

/// This test runs intersection with all portals of the toy detector
// TODO: use runge-kutta stepping
TEST(ALGEBRA_PLUGIN, geometry_discovery) {
    using namespace navigation;

    vecmem::host_memory_resource host_mr;
    auto det = create_toy_geometry(host_mr);

    // Create the navigator
    using inspector_t = aggregate_inspector<object_tracer<status::e_on_target>,
                                            print_inspector>;
    using navigator_t = navigator<decltype(det), inspector_t>;
    using stepper_t = line_stepper<free_track_parameters>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

    // Propagator
    propagator_t prop(stepper_t{}, navigator_t{det});

    unsigned int theta_steps = 1;
    unsigned int phi_steps = 1;

    const point3 ori{0., 0., 0.};
    // d.volume_by_pos(ori).index();

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
            ray r{ori, dir};
            const auto intersection_trace = shoot_ray(det, r);

            free_track_parameters track(ori, 0, dir, -1);
            propagator_t::state prop_state(track);

            prop.propagate(prop_state);

            // Retrieve navigation information
            auto &obj_tracer =
                prop_state._navigation.inspector()
                    .template get<object_tracer<status::e_on_target>>();
            auto &debug_printer = prop_state._navigation.inspector()
                                      .template get<print_inspector>();

            // Compare intersection records
            EXPECT_EQ(obj_tracer.object_trace.size(),
                      intersection_trace.size());

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
