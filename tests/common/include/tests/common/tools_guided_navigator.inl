/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/tools/line_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "detray/tools/rk_stepper.hpp"
#include "detray/tools/track.hpp"
#include "detray/tools/volume_graph.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"

/// @note __plugin has to be defined with a preprocessor command
namespace detray {

using vector3 = __plugin::vector3<detray::scalar>;

/** A navigation inspector that relays information about the encountered
 *  objects the way we need them to compare with the ray
 */
template <int navigation_status = 0,
          template <typename...> class vector_t = dvector>
struct object_tracer {

    // record all object id the navigator encounters
    vector_t<intersection> object_trace = {};

    template <typename state_type>
    auto operator()(state_type &state, const char * /*message*/) {
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
    auto operator()(state_type &state, const char *message) {
        std::string msg(message);

        debug_stream << msg << std::endl;

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
    auto operator()(state_type &state, const char *message) {
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

}  // namespace detray

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, guided_navigator) {
    using namespace detray;

    using inspector_t = aggregate_inspector<object_tracer<1>, print_inspector>;
    using b_field_t = constant_magnetic_field<>;
    using stepper_t = rk_stepper<b_field_t, free_track_parameters>;

    constexpr scalar tol = 1e-4;

    vecmem::host_memory_resource host_mr;

    point3 pos{0., 0., 0.};
    vector3 mom{1., 0., 0.};
    free_track_parameters pilot_track(pos, 0, mom, -1);

    vector3 B{0, 0, 2 * unit_constants::T};
    b_field_t b_field(B);

    stepper_t stepper(b_field);
    stepper_t::state step_state(pilot_track);

    // Number of plane surfaces
    dindex n_surfaces = 10;
    // Total distance between the all surfaces, as seen by the stepper
    scalar tel_length = 500. * unit_constants::mm;
    // Build telescope detector
    const auto telescope_det = create_telescope_detector(
        host_mr, pilot_track, stepper, n_surfaces, tel_length);

    // Building the telescope geometry should not change the track state
    ASSERT_TRUE(pilot_track.pos() == pos);
    // Build the graph
    volume_graph graph(telescope_det);
    std::cout << graph.to_string() << std::endl;

    using guided_navigator = navigator<decltype(telescope_det), inspector_t>;

    guided_navigator nav(telescope_det);
    guided_navigator::state nav_state;
    // Set initial volume (no grid yet)
    nav_state.set_volume(0u);

    //
    // Use navigator to step through telescope
    //
    bool heartbeat = true;
    while (heartbeat) {
        // (Re-)target
        heartbeat &= nav.target(nav_state, step_state);
        // Take the step
        heartbeat &= stepper.step(step_state, nav_state());
        // Enforce evaluation of only the next surface in the telescope
        nav_state.set_trust_level(
            guided_navigator::navigation_trust_level::e_high_trust);
        // And check the status
        heartbeat &= nav.status(nav_state, step_state);
    }

    // Check that navigator exited
    auto &debug_printer = nav_state.inspector().template get<print_inspector>();
    ASSERT_TRUE(nav_state.status() == -3) << debug_printer.to_string();

    // sequence of surfaces we expect to see
    std::vector<dindex> sf_sequence = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    auto &obj_tracer = nav_state.inspector().template get<object_tracer<1>>();
    // Check the surfaces that have been visited by the navigation
    ASSERT_TRUE(obj_tracer.object_trace.size() == sf_sequence.size());
    for (size_t i = 0; i < n_surfaces; ++i) {
        const auto &candidate = obj_tracer.object_trace[i];
        EXPECT_TRUE(candidate.index == sf_sequence[i]);
    }

    //
    // Directly step through detector
    //
    // Reset track
    free_track_parameters new_track(pos, 0, mom, -1);
    stepper_t::state new_step_state(new_track);
    scalar step_size = tel_length / n_surfaces;
    heartbeat = true;
    for (size_t i = 0; i < n_surfaces - 1; ++i) {
        // Take the step
        heartbeat &= stepper.step(new_step_state, step_size);
    }
    EXPECT_TRUE(heartbeat);
    EXPECT_NEAR(step_state._path_accumulated, new_step_state._path_accumulated,
                tol);
    EXPECT_NEAR(getter::norm(new_track.pos() - pilot_track.pos()), 0, tol);
    EXPECT_NEAR(getter::norm(new_track.dir() - pilot_track.dir()), 0, tol);
}
