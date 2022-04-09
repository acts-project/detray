/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/propagator/aborters.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/navigation_policies.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/propagator/track.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"
#include "tests/common/tools/inspectors.hpp"

/// @note __plugin has to be defined with a preprocessor command
namespace detray {

using vector3 = __plugin::vector3<detray::scalar>;

}

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, guided_navigator) {
    using namespace detray;
    using namespace navigation;

    vecmem::host_memory_resource host_mr;

    // Module positions along z-axis
    std::vector<scalar> positions = {0.,  10., 20., 30., 40., 50.,
                                     60., 70,  80,  90., 100.};
    // Build telescope detector with unbounded planes
    const auto telescope_det =
        create_telescope_detector<unbounded>(host_mr, positions);

    // Inspectors are optional, of course
    using inspector_t = aggregate_inspector<object_tracer<status::e_on_target>,
                                            print_inspector>;
    using b_field_t = constant_magnetic_field<>;
    using runge_kutta_stepper =
        rk_stepper<b_field_t, free_track_parameters, guided_navigation>;
    using guided_navigator = navigator<decltype(telescope_det), inspector_t>;
    using actor_chain_t =
        actor_chain<dtuple, pathlimit_aborter, guided_navigation>;
    using propagator_t =
        propagator<runge_kutta_stepper, guided_navigator, actor_chain_t>;

    // track must point into the direction of the telescope
    point3 pos{0., 0., 0.};
    vector3 mom{0., 0., 1.};
    free_track_parameters track(pos, 0, mom, -1);
    vector3 B{0, 0, 1 * unit_constants::T};
    b_field_t b_field(B);

    runge_kutta_stepper stepper(b_field);
    guided_navigator nav(telescope_det);

    // Actors
    pathlimit_aborter::state_type pathlimit{1. * unit_constants::m};
    guided_navigation::state_type policy{};
    actor_chain_t::state actor_states = std::tie(pathlimit, policy);

    // Propagator
    propagator_t p(std::move(stepper), std::move(nav));
    propagator_t::state guided_state(track, actor_states);

    // Propagate
    p.propagate(guided_state);

    auto &nav_state = guided_state._navigation;
    auto &debug_printer = nav_state.inspector().template get<print_inspector>();
    auto &obj_tracer = nav_state.inspector()
                           .template get<object_tracer<status::e_on_target>>();

    // Check that navigator exited
    ASSERT_TRUE(nav_state.is_complete()) << debug_printer.to_string();

    // sequence of surface ids we expect to see
    std::vector<dindex> sf_sequence = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    // Check the surfaces that have been visited by the navigation
    EXPECT_EQ(obj_tracer.object_trace.size(), sf_sequence.size());
    for (size_t i = 0; i < sf_sequence.size(); ++i) {
        auto &candidate = obj_tracer.object_trace[i];
        EXPECT_TRUE(candidate.index == sf_sequence[i]);
    }
}
