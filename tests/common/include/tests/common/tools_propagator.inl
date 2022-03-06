/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/transform_store.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/io/csv_io.hpp"
#include "detray/tools/line_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "detray/tools/propagator.hpp"
#include "detray/tools/rk_stepper.hpp"
#include "detray/tools/track.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/helix_gun.hpp"
#include "tests/common/tools/read_geometry.hpp"

constexpr scalar epsilon = 1e-4;

// This tests the basic functionality of the propagator
TEST(ALGEBRA_PLUGIN, propagator_line_stepper) {
    vecmem::host_memory_resource host_mr;

    using namespace detray;
    using namespace __plugin;

    // auto [d, name_map] = read_from_csv(tml_files, host_mr);
    auto d = create_toy_geometry(host_mr);

    // Create the navigator
    using detray_navigator = navigator<decltype(d)>;
    using detray_track = free_track_parameters;

    __plugin::point3<scalar> pos{0., 0., 0.};
    __plugin::vector3<scalar> mom{1., 1., 0.};
    detray_track traj(pos, 0, mom, -1);

    using detray_stepper = line_stepper<detray_track>;

    detray_stepper s;
    detray_navigator n(std::move(d));

    using detray_propagator = propagator<detray_stepper, detray_navigator>;
    detray_propagator p(std::move(s), std::move(n));

    void_propagator_inspector vi;

    detray_propagator::state state(traj);

    /*auto end = */  // p.propagate(state, vi);
}

struct helix_inspector {

    helix_inspector(const helix_gun& helix) : _helix(helix) {}

    template <typename navigator_state_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE void operator()(const navigator_state_t& /*navigation*/,
                                       const stepper_state_t& stepping) {
        auto pos = stepping().pos();
        auto true_pos = _helix(stepping._path_accumulated);

        auto relative_error = 1 / stepping._path_accumulated * (pos - true_pos);

        EXPECT_NEAR(getter::norm(relative_error), 0, epsilon);
    }

    helix_gun _helix;
};

struct print_inspector {

    template <typename navigator_state_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE void operator()(const navigator_state_t& navigation,
                                       const stepper_state_t& stepping) {
        std::stringstream stream;

        stream << std::left << std::setw(30);
        switch (static_cast<int>(navigation.status())) {
            case -3:
                stream << "status: on_target";
                break;
            case -2:
                stream << "status: abort";
                break;
            case -1:
                stream << "status: unknowm";
                break;
            case 0:
                stream << "status: towards_surface";
                break;
            case 1:
                stream << "status: on_surface";
                break;
        };

        if (navigation.volume() == dindex_invalid) {
            stream << "volume: " << std::setw(10) << "invalid";
        } else {
            stream << "volume: " << std::setw(10) << navigation.volume();
        }

        if (navigation.on_object() == dindex_invalid) {
            stream << "surface: " << std::setw(14) << "invalid";
        } else {
            stream << "surface: " << std::setw(14) << navigation.on_object();
        }

        stream << "step_size: " << std::setw(10) << stepping._step_size;

        // std::cout << stream.str() << std::endl;
    }
};

struct combined_inspector {
    helix_inspector _hi;
    print_inspector _pi;

    template <typename navigator_state_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE void operator()(navigator_state_t& navigation,
                                       const stepper_state_t& stepping) {
        _hi(navigation, stepping);
        _pi(navigation, stepping);
    }
};

/** A navigation inspector that prints information about the current navigation
 * state. Meant for debugging.
 */
struct navigator_print_inspector {

    // Debug output if an error in the trace is discovered
    std::stringstream debug_stream;

    template <typename state_type>
    auto operator()(state_type& state, const char* message) {
        std::string msg(message);

        debug_stream << msg << std::endl;

        debug_stream << "Volume\t\t\t\t\t\t" << state.volume() << std::endl;
        debug_stream << "surface kernel size\t\t" << state.candidates().size()
                     << std::endl;

        debug_stream << "Surface candidates: " << std::endl;
        for (const auto& sf_cand : state.candidates()) {
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
    auto operator()(state_type& state, const char* message) {
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

class PropagatorWithRkStepper
    : public ::testing::TestWithParam<__plugin::vector3<scalar>> {};

TEST_P(PropagatorWithRkStepper, propagator_rk_stepper) {

    using namespace detray;
    using namespace __plugin;
    using point3 = __plugin::point3<scalar>;
    using vector3 = __plugin::vector3<scalar>;

    // geomery navigation configurations
    constexpr unsigned int theta_steps = 100;
    constexpr unsigned int phi_steps = 100;

    // detector configuration
    constexpr std::size_t n_brl_layers = 4;
    constexpr std::size_t n_edc_layers = 7;

    vecmem::host_memory_resource host_mr;
    auto d = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    // Create the navigator
    using detray_navigator = navigator<decltype(d), navigator_print_inspector>;
    using detray_b_field = constant_magnetic_field<>;
    using detray_stepper = rk_stepper<detray_b_field, free_track_parameters>;
    using detray_propagator = propagator<detray_stepper, detray_navigator>;

    // Constant magnetic field
    // vector3 B{0, 0, 2 * unit_constants::T};
    vector3 B = GetParam();
    detray_b_field b_field(B);

    // Set initial position and momentum of tracks
    const point3 ori{0., 0., 0.};

    detray_stepper s(b_field);
    detray_navigator n(std::move(d));
    detray_propagator p(std::move(s), std::move(n));

    // Loops of theta values ]0,pi[
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.001 + itheta * (M_PI - 0.001) / theta_steps;
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        // Loops of phi values [-pi, pi]
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {
            // The direction
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            scalar sin_phi = std::sin(phi);
            scalar cos_phi = std::cos(phi);
            vector3 mom{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};
            mom = 10. * mom;

            // intialize a track
            free_track_parameters traj(ori, 0, mom, -1);
            helix_gun helix(traj, B);
            combined_inspector ci{helix_inspector(helix), print_inspector{}};

            detray_propagator::state state(traj);

            p.propagate(state, ci);

            // Ensure that the tracks reach the end of the world
            EXPECT_EQ(state._navigation.volume(), dindex_invalid);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(PropagatorValidation, PropagatorWithRkStepper,
                         ::testing::Values(__plugin::vector3<scalar>{
                             0, 0, 2 * unit_constants::T}));

/*
INSTANTIATE_TEST_SUITE_P(
    PropagatorValidation, PropagatorWithRkStepper,
    ::testing::Values(__plugin::vector3<scalar>{0, 0, 2 * unit_constants::T},
__plugin::vector3<scalar>{1 * unit_constants::T, 1 * unit_constants::T, 1 *
unit_constants::T}));
*/
