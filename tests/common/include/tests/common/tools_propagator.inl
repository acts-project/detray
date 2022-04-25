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
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/helix_gun.hpp"
#include "tests/common/tools/inspectors.hpp"

using namespace detray;

namespace {

constexpr scalar epsilon = 5e-4;
constexpr scalar path_limit = 5 * unit_constants::cm;

/// Compare helical track positions for stepper
struct helix_inspector : actor {

    /// Keeps the state of a helix gun to calculate track positions
    struct helix_inspector_state {
        helix_inspector_state(helix_gun &&h) : _helix(h) {}
        helix_gun _helix;
    };

    using size_type = __plugin::size_type;
    using matrix_operator = standard_matrix_operator<scalar>;
    using state_type = helix_inspector_state;

    /// Check that the stepper remains on the right helical track for its pos.
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        const state_type &inspector_state,
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

    using namespace detray;
    using namespace propagation;

    // geomery navigation configurations
    constexpr unsigned int theta_steps = 50;
    constexpr unsigned int phi_steps = 50;

    // detector configuration
    constexpr std::size_t n_brl_layers = 4;
    constexpr std::size_t n_edc_layers = 7;

    vecmem::host_memory_resource host_mr;
    auto d = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    // Create the navigator
    using navigator_t = navigator<decltype(d)>;
    using b_field_t = constant_magnetic_field<>;
    using constraints_t = constrained_step<>;
    using policy_t = stepper_default_policy;
    using stepper_t =
        rk_stepper<b_field_t, free_track_parameters, constraints_t, policy_t>;
    using actor_chain_t = actor_chain<dtuple, helix_inspector, print_inspector,
                                      pathlimit_aborter>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Constant magnetic field
    vector3 B = GetParam();
    b_field_t b_field(B);

    // Propagator is built from the stepper and navigator
    propagator_t p(stepper_t{b_field}, navigator_t{d});

    // Set origin position of tracks
    const point3 ori{0., 0., 0.};

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

            // intialize a track
            vector3 mom{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};
            vector::normalize(mom);
            mom = 10. * unit_constants::GeV * mom;
            free_track_parameters traj(ori, 0, mom, -1);
            free_track_parameters lim_traj(ori, 0, mom, -1);
            traj.set_overstep_tolerance(-10 * unit_constants::um);
            lim_traj.set_overstep_tolerance(-10 * unit_constants::um);

            // Build actor states: the helix inspector can be shared
            helix_inspector::state_type helix_insp_state{helix_gun{traj, &B}};
            print_inspector::state_type print_insp_state{};
            print_inspector::state_type lim_print_insp_state{};
            pathlimit_aborter::state_type unlimted_aborter_state{};
            pathlimit_aborter::state_type pathlimit_aborter_state{path_limit};

            // Create actor states tuples
            actor_chain_t::state actor_states = std::tie(
                helix_insp_state, print_insp_state, unlimted_aborter_state);
            actor_chain_t::state lim_actor_states =
                std::tie(helix_insp_state, lim_print_insp_state,
                         pathlimit_aborter_state);

            // Init propagator states
            propagator_t::state state(traj, actor_states);
            propagator_t::state lim_state(lim_traj, lim_actor_states);

            // Set step constraints
            state._stepping
                .template set_constraint<step::constraint::e_accuracy>(
                    5. * unit_constants::mm);
            lim_state._stepping
                .template set_constraint<step::constraint::e_accuracy>(
                    5. * unit_constants::mm);

            // Propagate the entire detector
            ASSERT_TRUE(p.propagate(state))
                << print_insp_state.to_string() << std::endl;

            // Propagate with path limit
            ASSERT_NEAR(pathlimit_aborter_state.path_limit(), path_limit,
                        epsilon);
            ASSERT_FALSE(p.propagate(lim_state))
                << lim_print_insp_state.to_string() << std::endl;
            ASSERT_TRUE(lim_state._stepping.path_length() <
                        path_limit + epsilon);
        }
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
