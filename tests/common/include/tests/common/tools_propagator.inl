/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/propagator/aborters.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/propagator/track.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/inspectors.hpp"
#include "tests/common/tools/test_trajectories.hpp"
#include "tests/common/tools/track_generators.hpp"

using namespace detray;
using matrix_operator = standard_matrix_operator<scalar>;
using mag_field_t = constant_magnetic_field<>;

namespace {

constexpr scalar epsilon = 5e-4;
constexpr scalar path_limit = 5 * unit_constants::cm;

/// Compare helical track positions for stepper
struct helix_inspector : actor {

    /// Keeps the state of a helix gun to calculate track positions
    struct state {
        state(helix &&h) : _helix(h) {}
        helix _helix;
    };

    using size_type = __plugin::size_type;
    using matrix_operator = standard_matrix_operator<scalar>;

    /// Check that the stepper remains on the right helical track for its pos.
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        const state &inspector_state,
        const propagator_state_t &prop_state) const {

        const auto &stepping = prop_state._stepping;
        const auto pos = stepping().pos();
        const auto true_pos = inspector_state._helix(stepping.path_length());

        const point3 relative_error{1 / stepping.path_length() *
                                    (pos - true_pos)};

        ASSERT_NEAR(getter::norm(relative_error), 0, epsilon);

        auto true_J = inspector_state._helix.jacobian(stepping.path_length());
        for (size_type i = 0; i < e_free_size; i++) {
            for (size_type j = 0; j < e_free_size; j++) {
                ASSERT_NEAR(
                    matrix_operator().element(stepping._jac_transport, i, j),
                    matrix_operator().element(true_J, i, j),
                    stepping.path_length() * epsilon * 10);
            }
        }
    }
};

}  // anonymous namespace

/// This tests the basic functionality of the propagator
TEST(ALGEBRA_PLUGIN, propagator_line_stepper) {

    vecmem::host_memory_resource host_mr;
    auto d = create_toy_geometry(host_mr);

    using navigator_t = navigator<decltype(d), navigation::print_inspector>;
    using track_t = free_track_parameters;
    using stepper_t = line_stepper<track_t>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

    point3 pos{0., 0., 0.};
    vector3 mom{1., 1., 0.};
    track_t traj(pos, 0, mom, -1);

    propagator_t p(stepper_t{}, navigator_t{d});

    propagator_t::state state(traj);

    EXPECT_TRUE(p.propagate(state))
        << state._navigation.inspector().to_string() << std::endl;
}

class PropagatorWithRkStepper
    : public ::testing::TestWithParam<__plugin::vector3<scalar>> {};

TEST_P(PropagatorWithRkStepper, propagator_rk_stepper) {

    // geomery navigation configurations
    constexpr unsigned int theta_steps = 50;
    constexpr unsigned int phi_steps = 50;

    // Set origin position of tracks
    const point3 ori{0., 0., 0.};
    constexpr scalar mom = 10. * unit_constants::GeV;

    // detector configuration
    constexpr std::size_t n_brl_layers = 4;
    constexpr std::size_t n_edc_layers = 7;

    vecmem::host_memory_resource host_mr;
    auto d = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    // Create the navigator
    using track_t = free_track_parameters;
    using navigator_t = navigator<decltype(d)>;
    using constraints_t = constrained_step<>;
    using stepper_t = rk_stepper<mag_field_t, track_t, constraints_t>;
    using actor_chain_t =
        actor_chain<dtuple, helix_inspector, propagation::print_inspector,
                    pathlimit_aborter>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Constant magnetic field
    vector3 B = GetParam();
    mag_field_t b_field(B);

    // Propagator is built from the stepper and navigator
    propagator_t p(stepper_t{b_field}, navigator_t{d});

    // Iterate through uniformly distributed momentum directions
    for (auto traj :
         uniform_track_generator<track_t>(theta_steps, phi_steps, ori, mom)) {
        // Genrate track state used for propagation with pathlimit
        free_track_parameters lim_traj(traj);

        traj.set_overstep_tolerance(-10 * unit_constants::um);
        lim_traj.set_overstep_tolerance(-10 * unit_constants::um);

        // Build actor states: the helix inspector can be shared
        helix_inspector::state helix_insp_state{helix{traj, &B}};
        propagation::print_inspector::state print_insp_state{};
        propagation::print_inspector::state lim_print_insp_state{};
        pathlimit_aborter::state unlimted_aborter_state{};
        pathlimit_aborter::state pathlimit_aborter_state{path_limit};

        // Create actor states tuples
        actor_chain_t::state actor_states = std::tie(
            helix_insp_state, print_insp_state, unlimted_aborter_state);
        actor_chain_t::state lim_actor_states = std::tie(
            helix_insp_state, lim_print_insp_state, pathlimit_aborter_state);

        // Init propagator states
        propagator_t::state state(traj, actor_states);
        propagator_t::state lim_state(lim_traj, lim_actor_states);

        // Set step constraints
        state._stepping.template set_constraint<step::constraint::e_accuracy>(
            5. * unit_constants::mm);
        lim_state._stepping
            .template set_constraint<step::constraint::e_accuracy>(
                5. * unit_constants::mm);

        // Propagate the entire detector
        ASSERT_TRUE(p.propagate(state))
            << print_insp_state.to_string() << std::endl;

        // Propagate with path limit
        ASSERT_NEAR(pathlimit_aborter_state.path_limit(), path_limit, epsilon);
        ASSERT_FALSE(p.propagate(lim_state))
            << lim_print_insp_state.to_string() << std::endl;
        ASSERT_TRUE(lim_state._stepping.path_length() < path_limit + epsilon);
    }
}

INSTANTIATE_TEST_SUITE_P(PropagatorValidation1, PropagatorWithRkStepper,
                         ::testing::Values(__plugin::vector3<scalar>{
                             0. * unit_constants::T, 0. * unit_constants::T,
                             2. * unit_constants::T}));

INSTANTIATE_TEST_SUITE_P(PropagatorValidation2, PropagatorWithRkStepper,
                         ::testing::Values(__plugin::vector3<scalar>{
                             0. * unit_constants::T, 1. * unit_constants::T,
                             1. * unit_constants::T}));

INSTANTIATE_TEST_SUITE_P(PropagatorValidation3, PropagatorWithRkStepper,
                         ::testing::Values(__plugin::vector3<scalar>{
                             1. * unit_constants::T, 0. * unit_constants::T,
                             1. * unit_constants::T}));

INSTANTIATE_TEST_SUITE_P(PropagatorValidation4, PropagatorWithRkStepper,
                         ::testing::Values(__plugin::vector3<scalar>{
                             1. * unit_constants::T, 1. * unit_constants::T,
                             1. * unit_constants::T}));