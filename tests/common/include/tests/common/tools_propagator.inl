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
#include "detray/intersection/detail/trajectories.hpp"
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
#include "tests/common/tools/track_generators.hpp"

using namespace detray;

namespace {

constexpr scalar epsilon{1e-3};
constexpr scalar path_limit{5 * unit_constants::cm};

/// Compare helical track positions for stepper
struct helix_inspector : actor {

    /// Keeps the state of a helix gun to calculate track positions
    struct state {
        state(detail::helix &&h) : _helix(h) {}
        detail::helix _helix;
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

        // Nothing has happened yet (first call of actor chain)
        if (stepping.path_length() < epsilon) {
            return;
        }
        const point3 relative_error{scalar{1.} / stepping.path_length() *
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

/// Test basic functionality of the propagator using a straight line stepper
TEST(ALGEBRA_PLUGIN, propagator_line_stepper) {

    vecmem::host_memory_resource host_mr;
    const auto d = create_toy_geometry(host_mr);

    using navigator_t = navigator<decltype(d), navigation::print_inspector>;
    using stepper_t = line_stepper<free_track_parameters>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

    const point3 pos{0., 0., 0.};
    const vector3 mom{1., 1., 0.};
    free_track_parameters track(pos, 0, mom, -1);

    propagator_t p(stepper_t{}, navigator_t{d});

    propagator_t::state state(track);

    EXPECT_TRUE(p.propagate(state))
        << state._navigation.inspector().to_string() << std::endl;
}

class PropagatorWithRkStepper
    : public ::testing::TestWithParam<
          std::tuple<__plugin::vector3<scalar>, scalar, scalar>> {};

/// Test propagation in a magnetic field using a Runge-Kutta stepper
TEST_P(PropagatorWithRkStepper, propagator_rk_stepper) {

    // geomery navigation configurations
    constexpr unsigned int theta_steps{50};
    constexpr unsigned int phi_steps{50};

    // Set origin position of tracks
    const point3 ori{0., 0., 0.};
    constexpr scalar mom{10. * unit_constants::GeV};

    // detector configuration
    constexpr std::size_t n_brl_layers{4};
    constexpr std::size_t n_edc_layers{7};
    vecmem::host_memory_resource host_mr;
    const auto d = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    // Create the navigator
    // using navigator_t = navigator<decltype(d), navigation::print_inspector>;
    using navigator_t = navigator<decltype(d)>;
    using b_field_t = constant_magnetic_field<>;
    using track_t = free_track_parameters;
    using constraints_t = constrained_step<>;
    using policy_t = stepper_default_policy;
    using stepper_t = rk_stepper<b_field_t, track_t, constraints_t, policy_t>;
    using actor_chain_t =
        actor_chain<dtuple, helix_inspector, propagation::print_inspector,
                    pathlimit_aborter>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Test parameters
    const vector3 B = std::get<0>(GetParam());
    const scalar overstep_tol = std::get<1>(GetParam());
    const scalar step_constr = std::get<2>(GetParam());
    // Constant magnetic field
    const b_field_t b_field(B);

    // Propagator is built from the stepper and navigator
    propagator_t p(stepper_t{b_field}, navigator_t{d});

    // Iterate through uniformly distributed momentum directions
    for (auto track :
         uniform_track_generator<track_t>(theta_steps, phi_steps, ori, mom)) {
        // Genrate second track state used for propagation with pathlimit
        track_t lim_track(track);

        track.set_overstep_tolerance(overstep_tol);
        lim_track.set_overstep_tolerance(overstep_tol);

        // Build actor states: the helix inspector can be shared
        helix_inspector::state helix_insp_state{detail::helix{track, &B}};
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
        propagator_t::state state(track, actor_states);
        propagator_t::state lim_state(lim_track, lim_actor_states);

        // Set step constraints
        state._stepping.template set_constraint<step::constraint::e_accuracy>(
            step_constr);
        lim_state._stepping
            .template set_constraint<step::constraint::e_accuracy>(step_constr);

        // Propagate the entire detector
        ASSERT_TRUE(p.propagate(state))
            << print_insp_state.to_string() << std::endl;
        //<< state._navigation.inspector().to_string() << std::endl;

        // Propagate with path limit
        ASSERT_NEAR(pathlimit_aborter_state.path_limit(), path_limit, epsilon);
        ASSERT_FALSE(p.propagate(lim_state))
            << lim_print_insp_state.to_string() << std::endl;
        //<< lim_state._navigation.inspector().to_string() << std::endl;

        ASSERT_TRUE(lim_state._stepping.path_length() <
                    std::abs(path_limit) + epsilon)
            << "path length: " << lim_state._stepping.path_length()
            << ", path limit: " << path_limit << std::endl;
        //<< state._navigation.inspector().to_string() << std::endl;
    }
}

// Realistic case
INSTANTIATE_TEST_SUITE_P(PropagatorValidation1, PropagatorWithRkStepper,
                         ::testing::Values(std::make_tuple(
                             __plugin::vector3<scalar>{0. * unit_constants::T,
                                                       0. * unit_constants::T,
                                                       2. * unit_constants::T},
                             -7. * unit_constants::um,
                             std::numeric_limits<scalar>::max())));

// Add some restrictions for more frequent navigation updates in the cases of
// non-z-aligned B-fields
INSTANTIATE_TEST_SUITE_P(PropagatorValidation2, PropagatorWithRkStepper,
                         ::testing::Values(std::make_tuple(
                             __plugin::vector3<scalar>{0. * unit_constants::T,
                                                       1. * unit_constants::T,
                                                       1. * unit_constants::T},
                             -10. * unit_constants::um,
                             5. * unit_constants::mm)));

INSTANTIATE_TEST_SUITE_P(PropagatorValidation3, PropagatorWithRkStepper,
                         ::testing::Values(std::make_tuple(
                             __plugin::vector3<scalar>{1. * unit_constants::T,
                                                       0. * unit_constants::T,
                                                       1. * unit_constants::T},
                             -10. * unit_constants::um,
                             5. * unit_constants::mm)));

INSTANTIATE_TEST_SUITE_P(PropagatorValidation4, PropagatorWithRkStepper,
                         ::testing::Values(std::make_tuple(
                             __plugin::vector3<scalar>{1. * unit_constants::T,
                                                       1. * unit_constants::T,
                                                       1. * unit_constants::T},
                             -10. * unit_constants::um,
                             5. * unit_constants::mm)));
