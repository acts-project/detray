/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <sstream>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/field/constant_magnetic_field.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/inspectors.hpp"
#include "tests/common/tools/particle_gun.hpp"
#include "tests/common/tools/track_generators.hpp"

using namespace detray;

/// This test runs intersection with all portals of the toy detector with a ray
/// and then compares the intersection trace with a straight line navigation.
TEST(ALGEBRA_PLUGIN, straight_line_navigation) {
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

        // Shoot ray through the detector and record all surfaces it encounters
        const auto intersection_trace =
            particle_gun::shoot_particle(det, test_ray);

        // Now follow that ray with a track and check, if we find the same
        // volumes and distances along the way
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
        for (std::size_t i = 0; i < intersection_trace.size(); ++i) {
            if (obj_tracer[i].index != intersection_trace[i].first) {
                // Intersection record at portal bound might be flipped
                // (the portals overlap completely)
                if (obj_tracer[i].index == intersection_trace[i + 1].first and
                    obj_tracer[i + 1].index == intersection_trace[i].first) {
                    // Have already checked the next record
                    ++i;
                    continue;
                }
            }
            EXPECT_EQ(obj_tracer[i].index, intersection_trace[i].first)
                << debug_printer.to_string() << debug_stream.str();
        }
    }
}

TEST(ALGEBRA_PLUGIN, helix_navigation) {
    using namespace navigation;

    vecmem::host_memory_resource host_mr;
    auto det = create_toy_geometry(host_mr);

    // Create the navigator
    using inspector_t = aggregate_inspector<object_tracer<status::e_on_target>,
                                            print_inspector>;
    using navigator_t = navigator<decltype(det), inspector_t>;
    using mag_field_t = constant_magnetic_field<>;
    using stepper_t = rk_stepper<mag_field_t, free_track_parameters>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    2. * unit_constants::T};
    mag_field_t b_field(B);

    // Propagator
    propagator_t prop(stepper_t{b_field}, navigator_t{det});

    unsigned int theta_steps = 1;
    unsigned int phi_steps = 1;

    const point3 ori{0., 0., 0.};

    // Iterate through uniformly distributed momentum directions
    for (const auto test_track : uniform_track_generator<free_track_parameters>(
             theta_steps, phi_steps, ori)) {
        helix h2(test_track, &B);
    }
}