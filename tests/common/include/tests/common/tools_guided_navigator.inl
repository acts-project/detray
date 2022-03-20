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

    using inspector_t = aggregate_inspector<object_tracer<1>, print_inspector>;
    using b_field_t = constant_magnetic_field<>;
    using stepper_t = rk_stepper<b_field_t, free_track_parameters>;

    constexpr scalar tol = 1e-4;

    vecmem::host_memory_resource host_mr;

    point3 pos{0., 0., 0.};
    vector3 mom{1., 0., 0.};
    free_track_parameters pilot_track(pos, 0, mom, -1);
    pilot_track.set_overstep_tolerance(-10 * unit_constants::um);

    vector3 B{0, 0, 2 * unit_constants::T};
    b_field_t b_field(B);

    stepper_t stepper(b_field);
    stepper_t::state step_state(pilot_track);

    // Create telescope detector of either unbounded planes or rectangles
    constexpr bool unbounded = true;
    // constexpr bool rectangles = false;
    //  Number of plane surfaces
    dindex n_surfaces = 10;
    // Total distance between all of the surfaces, as seen by the stepper
    scalar tel_length = 500. * unit_constants::mm;

    std::vector<scalar> positions = {0.,   50., 100., 150., 200., 250.,
                                     300., 350, 400,  450., 500.};
    // Build telescope detector with unbounded planes
    const auto telescope_det = create_telescope_detector<unbounded>(
        host_mr, n_surfaces, tel_length, {}, pilot_track, stepper);
    // Build telescope detector with rectangular planes
    // const auto telescope_det = create_telescope_detector<rectangles>(
    //     host_mr, n_surfaces, tel_length, {}, pilot_track, stepper);

    using guided_navigator = navigator<decltype(telescope_det), inspector_t>;

    guided_navigator nav(telescope_det);
    guided_navigator::state nav_state;

    //
    // Use navigator to step through telescope
    //
    bool heartbeat = nav.init(nav_state, step_state);

    while (heartbeat) {
        // Take the step
        heartbeat &= stepper.step(step_state, nav_state);
        // Enforce evaluation of only the next surface in the telescope
        nav_state.set_trust_level(e_high_trust);
        // And check the status
        heartbeat &= nav.update(nav_state, step_state);
    }

    // Check that navigator exited
    auto &debug_printer = nav_state.inspector().template get<print_inspector>();
    ASSERT_TRUE(nav_state.is_complete()) << debug_printer.to_string();

    // sequence of surfaces we expect to see
    std::vector<dindex> sf_sequence = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    auto &obj_tracer = nav_state.inspector().template get<object_tracer<1>>();
    // Check the surfaces that have been visited by the navigation
    EXPECT_EQ(obj_tracer.object_trace.size(), sf_sequence.size());
    for (size_t i = 0; i < n_surfaces; ++i) {
        const auto &candidate = obj_tracer.object_trace[i];
        std::cout << candidate.index << std::endl;
        EXPECT_TRUE(candidate.index == sf_sequence[i]);
    }
    std::cout << obj_tracer.object_trace.back().index << std::endl;

    //
    // Step through detector directly
    //
    // dummy navigation struct
    struct navigation_state {
        scalar operator()() const { return _step_size; }
        inline void set_full_trust() {}
        inline void set_high_trust() {}
        inline void set_fair_trust() {}
        inline void set_no_trust() {}
        inline bool abort() { return false; }

        scalar _step_size;
    };
    navigation_state new_nav_state{tel_length / (n_surfaces - 1)};

    // Reset track
    free_track_parameters new_track(pos, 0, mom, -1);
    stepper_t::state new_step_state(new_track);

    heartbeat = true;
    for (size_t i = 0; i < n_surfaces - 1; ++i) {
        // Take a step
        heartbeat &= stepper.step(new_step_state, new_nav_state);
    }
    EXPECT_TRUE(heartbeat);
    EXPECT_NEAR(
        std::fabs(step_state._path_length - new_step_state._path_length) /
            step_state._path_length,
        0., tol);
    EXPECT_NEAR(getter::norm(new_track.pos() - pilot_track.pos()) /
                    getter::norm(pilot_track.pos()),
                0., tol);
    EXPECT_NEAR(getter::norm(new_track.dir() - pilot_track.dir()) /
                    getter::norm(pilot_track.dir()),
                0., tol);
}
