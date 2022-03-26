/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/tools/line_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "detray/tools/rk_stepper.hpp"
#include "detray/tools/track.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"
#include "tests/common/tools/inspectors.hpp"

/// @note __plugin has to be defined with a preprocessor command
namespace detray {

using vector3 = __plugin::vector3<detray::scalar>;
}

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, guided_navigator) {
    using namespace detray;
    using namespace detray::navigation;

    using inspector_t = aggregate_inspector<object_tracer<status::e_on_target>,
                                            print_inspector>;
    using b_field_t = constant_magnetic_field<>;
    using stepper_t = rk_stepper<b_field_t, free_track_parameters>;

    vecmem::host_memory_resource host_mr;

    //  Number of plane surfaces
    dindex n_surfaces = 11;
    // Module positions along z-axis
    std::vector<scalar> positions = {0.,  10., 20., 30., 40., 50.,
                                     60., 70,  80,  90., 100.};
    // Build telescope detector with unbounded planes
    const auto telescope_det =
        create_telescope_detector<unbounded>(host_mr, positions);

    using guided_navigator = navigator<decltype(telescope_det), inspector_t>;
    guided_navigator nav(telescope_det);
    guided_navigator::state nav_state;

    // track must point into the direction of the telescope
    point3 pos{0., 0., 0.};
    vector3 mom{0., 0., 1.};
    free_track_parameters track(pos, 0, mom, -1);

    vector3 B{0, 0, 1 * unit_constants::T};
    b_field_t b_field(B);

    stepper_t stepper(b_field);
    stepper_t::state step_state(track);

    //
    // Use navigator to step through telescope (re-evaluate only next surface)
    //
    bool heartbeat = nav.init(nav_state, step_state);

    while (heartbeat) {

        heartbeat &= stepper.step(step_state, nav_state);

        // Enforce evaluation of only the next surface in the telescope
        nav_state.set_trust_level(trust_level::e_high);

        heartbeat &= nav.update(nav_state, step_state);
    }

    auto &debug_printer = nav_state.inspector().template get<print_inspector>();
    auto &obj_tracer = nav_state.inspector()
                           .template get<object_tracer<status::e_on_target>>();

    // Check that navigator exited
    ASSERT_TRUE(nav_state.is_complete()) << debug_printer.to_string();

    // sequence of surfaces we expect to see
    std::vector<dindex> sf_sequence = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    // Check the surfaces that have been visited by the navigation
    EXPECT_EQ(obj_tracer.object_trace.size(), sf_sequence.size());
    for (size_t i = 0; i < n_surfaces; ++i) {
        const auto &candidate = obj_tracer.object_trace[i];
        EXPECT_TRUE(candidate.index == sf_sequence[i]);
    }
}
