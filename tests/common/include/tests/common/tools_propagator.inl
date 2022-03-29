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
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/propagator/track.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/helix_gun.hpp"
#include "tests/common/tools/inspectors.hpp"
#include "tests/common/tools/read_geometry.hpp"

constexpr scalar epsilon = 5e-4;
constexpr scalar path_limit = 2 * unit_constants::m;

struct print_actor : actor<> {

    using actor_type = print_actor;

    /// Printer id
    struct state {
        int printer_id = 1;
    };

    /// Do this actors work and then call the observers on the updated status
    void operator()(result<actor_type> &res) {
        std::cout << "This is print actor: " << res.state.printer_id
                  << std::endl;
    }

    template <typename subj_result_t>
    void operator()(result<actor_type> &res, subj_result_t &subject_result) {
        std::cout << "This is print actor: " << res.state.printer_id
                  << ", observing actor: " << subject_result.state.printer_id
                  << std::endl;
    }
};

template <typename... observers>
struct example_actor : actor<observers...> {

    using actor_type = example_actor<observers...>;

    /// Printer id
    struct state {
        int printer_id = 3;
    };

    /// Print your id
    void operator()(result<actor_type> &res) {
        std::cout << "This is something actor: " << res.state.printer_id
                  << std::endl;
    }

    /// Print your id
    template <typename result_t>
    void operator()(result<actor_type> &res, result_t & /*subject_result*/) {
        std::cout << "This is something actor: " << res.state.printer_id
                  << std::endl;
    }
};

// Test the actor chain on some dummy actor types
TEST(ALGEBRA_PLUGIN, actor_chain) {
    using namespace detray;

    // observer types
    using composite1 = composite_actor<dtuple, example_actor, print_actor>;
    using composite2 =
        composite_actor<dtuple, example_actor, print_actor, print_actor>;
    // using composite3 = composite_actor<do_smth_actor, dtuple, composite2>;

    using actor_chain_t = actor_chain<dtuple, composite1, composite2>;

    actor_chain_t actors{};
    actor_chain_t::state actor_states;

    actors(actor_states);
}

// This tests the basic functionality of the propagator
TEST(ALGEBRA_PLUGIN, propagator_line_stepper) {
    vecmem::host_memory_resource host_mr;

    using namespace detray;
    using namespace __plugin;

    // auto [d, name_map] = read_from_csv(tml_files, host_mr);
    auto d = create_toy_geometry(host_mr);

    // Create the navigator
    using navigator_t = navigator<decltype(d), navigation::print_inspector>;
    using track_t = free_track_parameters;

    __plugin::point3<scalar> pos{0., 0., 0.};
    __plugin::vector3<scalar> mom{1., 1., 0.};
    track_t traj(pos, 0, mom, -1);

    using stepper_t = line_stepper<track_t>;

    stepper_t s;
    navigator_t n(d);

    using propagator_t = propagator<stepper_t, navigator_t>;
    propagator_t p(std::move(s), std::move(n));

    propagation::void_inspector vi;

    propagator_t::state state(traj);
    state._stepping.set_path_limit(path_limit);

    EXPECT_TRUE(p.propagate(state, vi))
        << state._navigation.inspector().to_string() << std::endl;
}

struct helix_inspector {

    helix_inspector(const helix_gun &helix) : _helix(helix) {}

    template <typename navigator_state_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE void operator()(const navigator_state_t & /*navigation*/,
                                       const stepper_state_t &stepping) {
        auto pos = stepping().pos();
        const scalar path_accumulated =
            path_limit - stepping.dist_to_path_limit();
        auto true_pos = _helix(path_accumulated);

        __plugin::vector3<scalar> relative_error =
            1 / path_accumulated * (pos - true_pos);

        ASSERT_NEAR(getter::norm(relative_error), 0, epsilon);
    }

    helix_gun _helix;
};

struct combined_inspector {
    helix_inspector _hi;
    propagation::print_inspector _pi;

    template <typename navigator_state_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE void operator()(navigator_state_t &navigation,
                                       const stepper_state_t &stepping) {
        _hi(navigation, stepping);
        _pi(navigation, stepping);
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
    using navigator_t = navigator<decltype(d)>;
    using b_field_t = constant_magnetic_field<>;
    using constraints_t = constrained_step<>;
    using stepper_t =
        rk_stepper<b_field_t, free_track_parameters, constraints_t>;
    using propagator_t = propagator<stepper_t, navigator_t>;

    // Constant magnetic field
    vector3 B = GetParam();
    b_field_t b_field(B);

    stepper_t s(b_field);
    navigator_t n(d);
    propagator_t p(std::move(s), std::move(n));

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
            traj.set_overstep_tolerance(-10 * unit_constants::um);

            helix_gun helix(traj, B);
            combined_inspector ci{helix_inspector(helix),
                                  propagation::print_inspector{}};

            propagator_t::state state(traj);
            state._stepping.set_path_limit(path_limit);
            state._stepping
                .template set_constraint<step::constraint::e_accuracy>(
                    5. * unit_constants::mm);

            ASSERT_TRUE(p.propagate(state, ci))
                << ci._pi.to_string() << std::endl;
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
