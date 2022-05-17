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

    // Iterate through uniformly distributed momentum directions
    for (const auto test_ray :
         uniform_track_generator<ray>(theta_steps, phi_steps, ori)) {

        // Now follow that ray and check, if we find the same
        // volumes and distances along the way
        const auto intersection_trace = shoot_ray(det, test_ray);

        free_track_parameters track(test_ray.pos(), 0, test_ray.dir(), -1);
        propagator_t::state propagation(track);

        prop.propagate(propagation);

        // Retrieve navigation information
        auto &inspector = propagation._navigation.inspector();
        auto &obj_tracer =
            inspector.template get<object_tracer<status::e_on_target>>();
        auto &debug_printer = inspector.template get<print_inspector>();

        // Compare intersection records
        EXPECT_EQ(obj_tracer.object_trace.size(), intersection_trace.size());

        std::stringstream debug_stream;
        for (std::size_t intr_idx = 0; intr_idx < intersection_trace.size();
             ++intr_idx) {
            debug_stream << "-------Intersection trace\n"
                         << "ray gun: "
                         << "\tsf id: " << intersection_trace[intr_idx].first
                         << ", "
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
