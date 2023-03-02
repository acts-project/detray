/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/inspectors.hpp"

using namespace detray;
using transform3 = __plugin::transform3<detray::scalar>;

namespace {

constexpr scalar tol{1e-3f};
constexpr scalar path_limit{5.f * unit<scalar>::cm};

/// Compare helical track positions for stepper
struct helix_inspector : actor {

    /// Keeps the state of a helix gun to calculate track positions
    struct state {
        // navigation status for every step
        std::vector<navigation::status> _nav_status;
    };

    using matrix_operator = standard_matrix_operator<scalar>;
    using size_type = typename matrix_operator::size_ty;
    using scalar_type = typename matrix_operator::scalar_type;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using free_vector = matrix_type<e_free_size, 1>;

    // Kernel to get a free position at the last surface
    struct kernel {
        using output_type = free_vector;

        template <typename mask_group_t, typename index_t,
                  typename transform_store_t, typename surface_t,
                  typename stepper_state_t>
        DETRAY_HOST_DEVICE inline output_type operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform_store_t& trf_store, const surface_t& surface,
            const stepper_state_t& stepping) {

            const auto& trf3 = trf_store[surface.transform()];

            const auto& mask = mask_group[index];

            auto local_coordinate = mask.local_frame();

            return local_coordinate.bound_to_free_vector(
                trf3, mask, stepping._bound_params.vector());
        }
    };

    /// Check that the stepper remains on the right helical track for its pos.
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        state& inspector_state, const propagator_state_t& prop_state) const {

        const auto& navigation = prop_state._navigation;
        const auto& stepping = prop_state._stepping;

        inspector_state._nav_status.push_back(navigation.status());

        if (prop_state.param_type() == parameter_type::e_free) {
            return;
        }

        // Nothing has happened yet (first call of actor chain)
        if (stepping.path_length() < tol || stepping._s < tol) {
            return;
        }

        const auto& det = navigation.detector();
        const auto& surface_container = det->surfaces();
        const auto& trf_store = det->transform_store();
        const auto& mask_store = det->mask_store();

        // Surface
        const auto& surface =
            surface_container[stepping._bound_params.surface_link()];

        const auto free_vec = mask_store.template call<kernel>(
            surface.mask(), trf_store, surface, stepping);

        const auto last_pos =
            detail::track_helper<matrix_operator>().pos(free_vec);

        free_track_parameters<transform3> free_params;
        free_params.set_vector(free_vec);

        const auto bvec =
            stepping._magnetic_field.at(last_pos[0], last_pos[1], last_pos[2]);
        const vector3 b{bvec[0], bvec[1], bvec[2]};
        // const auto b = stepping._magnetic_field->get_field(last_pos, {});
        detail::helix<transform3> hlx(free_params, &b);

        const auto true_pos = hlx(stepping._s);

        const point3 relative_error{1.f / stepping._s *
                                    (stepping().pos() - true_pos)};

        ASSERT_NEAR(getter::norm(relative_error), 0.f, tol);

        auto true_J = hlx.jacobian(stepping._s);

        for (unsigned int i = 0u; i < e_free_size; i++) {
            for (unsigned int j = 0u; j < e_free_size; j++) {
                ASSERT_NEAR(
                    matrix_operator().element(stepping._jac_transport, i, j),
                    matrix_operator().element(true_J, i, j),
                    stepping._s * tol * 10.f);
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
    using stepper_t = line_stepper<transform3>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

    const point3 pos{0.f, 0.f, 0.f};
    const vector3 mom{1.f, 1.f, 0.f};
    free_track_parameters<transform3> track(pos, 0.f, mom, -1.f);

    propagator_t p(stepper_t{}, navigator_t{});

    propagator_t::state state(track, d);

    EXPECT_TRUE(p.propagate(state))
        << state._navigation.inspector().to_string() << std::endl;
}

class PropagatorWithRkStepper
    : public ::testing::TestWithParam<
          std::tuple<__plugin::vector3<scalar>, scalar, scalar>> {};

/// Test propagation in a magnetic field using a Runge-Kutta stepper
TEST_P(PropagatorWithRkStepper, propagator_rk_stepper) {

    // geomery navigation configurations
    constexpr std::size_t theta_steps{50u};
    constexpr std::size_t phi_steps{50u};

    // Set origin position of tracks
    const point3 ori{0.f, 0.f, 0.f};
    constexpr scalar mom{10.f * unit<scalar>::GeV};

    // detector configuration
    constexpr std::size_t n_brl_layers{4u};
    constexpr std::size_t n_edc_layers{7u};
    vecmem::host_memory_resource host_mr;

    // Construct the constant magnetic field.
    using b_field_t = decltype(
        create_toy_geometry(host_mr, n_brl_layers, n_edc_layers))::bfield_type;
    const vector3 B = std::get<0>(GetParam());

    const auto d = create_toy_geometry(
        host_mr,
        n_brl_layers, n_edc_layers);

    // Create the navigator
    // using navigator_t = navigator<decltype(d), navigation::print_inspector>;
    using navigator_t = navigator<decltype(d)>;
    using track_t = free_track_parameters<transform3>;
    using constraints_t = constrained_step<>;
    using policy_t = stepper_default_policy;
    using stepper_t =
        rk_stepper<b_field_t::view_t, transform3, constraints_t, policy_t>;
    using actor_chain_t =
        actor_chain<dtuple, helix_inspector, propagation::print_inspector,
                    pathlimit_aborter, parameter_transporter<transform3>,
                    pointwise_material_interactor<transform3>,
                    parameter_resetter<transform3>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Test parameters
    const scalar overstep_tol = std::get<1>(GetParam());
    const scalar step_constr = std::get<2>(GetParam());

    // Propagator is built from the stepper and navigator
    propagator_t p(stepper_t{}, navigator_t{});

    // Iterate through uniformly distributed momentum directions
    for (auto track :
         uniform_track_generator<track_t>(theta_steps, phi_steps, ori, mom)) {
        // Genrate second track state used for propagation with pathlimit
        track_t lim_track(track);

        track.set_overstep_tolerance(overstep_tol);
        lim_track.set_overstep_tolerance(overstep_tol);

        // Build actor states: the helix inspector can be shared
        helix_inspector::state helix_insp_state{};
        propagation::print_inspector::state print_insp_state{};
        propagation::print_inspector::state lim_print_insp_state{};
        pathlimit_aborter::state unlimted_aborter_state{};
        pathlimit_aborter::state pathlimit_aborter_state{path_limit};
        parameter_transporter<transform3>::state transporter_state{};
        pointwise_material_interactor<transform3>::state interactor_state{};
        parameter_resetter<transform3>::state resetter_state{};

        // Create actor states tuples
        auto actor_states =
            std::tie(helix_insp_state, print_insp_state, unlimted_aborter_state,
                     transporter_state, interactor_state, resetter_state);
        auto sync_actor_states =
            std::tie(helix_insp_state, print_insp_state, unlimted_aborter_state,
                     transporter_state, interactor_state, resetter_state);
        auto lim_actor_states = std::tie(
            helix_insp_state, lim_print_insp_state, pathlimit_aborter_state,
            transporter_state, interactor_state, resetter_state);

        // Init propagator states
        propagator_t::state state(track, d.get_bfield(), d);
        propagator_t::state lim_state(lim_track, d.get_bfield(), d);

        // Set step constraints
        state._stepping.template set_constraint<step::constraint::e_accuracy>(
            step_constr);
        lim_state._stepping
            .template set_constraint<step::constraint::e_accuracy>(step_constr);

        // Propagate the entire detector
        ASSERT_TRUE(p.propagate(state, actor_states))
            << print_insp_state.to_string() << std::endl;
        //<< state._navigation.inspector().to_string() << std::endl;

        // Propagate with path limit
        ASSERT_NEAR(pathlimit_aborter_state.path_limit(), path_limit, tol);
        ASSERT_FALSE(p.propagate(lim_state, lim_actor_states))
            << lim_print_insp_state.to_string() << std::endl;
        //<< lim_state._navigation.inspector().to_string() << std::endl;

        ASSERT_TRUE(lim_state._stepping.path_length() <
                    std::abs(path_limit) + tol)
            << "path length: " << lim_state._stepping.path_length()
            << ", path limit: " << path_limit << std::endl;
        //<< state._navigation.inspector().to_string() << std::endl;

        // Compare the navigation status vector between propagate and
        // propagate_sync function
        const auto nav_status =
            std::get<helix_inspector::state&>(actor_states)._nav_status;
        const auto sync_nav_status =
            std::get<helix_inspector::state&>(sync_actor_states)._nav_status;
        ASSERT_TRUE(nav_status.size() > 0);
        ASSERT_TRUE(nav_status == sync_nav_status);
    }
}

// Realistic case
INSTANTIATE_TEST_SUITE_P(
    PropagatorValidation1, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(
        __plugin::vector3<scalar>{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                  2.f * unit<scalar>::T},
        -7.f * unit<scalar>::um, std::numeric_limits<scalar>::max())));

// Add some restrictions for more frequent navigation updates in the cases of
// non-z-aligned B-fields
INSTANTIATE_TEST_SUITE_P(
    PropagatorValidation2, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(
        __plugin::vector3<scalar>{0.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                  1.f * unit<scalar>::T},
        -10.f * unit<scalar>::um, 5.f * unit<scalar>::mm)));

INSTANTIATE_TEST_SUITE_P(
    PropagatorValidation3, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(
        __plugin::vector3<scalar>{1.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                  1.f * unit<scalar>::T},
        -10.f * unit<scalar>::um, 5.f * unit<scalar>::mm)));

INSTANTIATE_TEST_SUITE_P(
    PropagatorValidation4, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(
        __plugin::vector3<scalar>{1.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                  1.f * unit<scalar>::T},
        -10.f * unit<scalar>::um, 5.f * unit<scalar>::mm)));
