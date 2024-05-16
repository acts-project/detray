/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/propagator/propagator.hpp"

#include "detray/definitions/units.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/detectors/build_toy_detector.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/navigation/detail/trajectories.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/inspectors.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using algebra_t = test::algebra;
using scalar_t = test::scalar;
using point3 = test::point3;
using vector3 = test::vector3;

namespace {

constexpr scalar_t tol{1e-3f};
constexpr scalar_t path_limit{5.f * unit<scalar_t>::cm};

/// Compare helical track positions for stepper
struct helix_inspector : actor {

    /// Keeps the state of a helix gun to calculate track positions
    struct state {
        /// navigation status for every step
        std::vector<navigation::status> _nav_status;
    };

    using matrix_operator = test::matrix_operator;

    /// Check that the stepper remains on the right helical track for its pos.
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        state& inspector_state, const propagator_state_t& prop_state) const {

        const auto& navigation = prop_state._navigation;
        const auto& stepping = prop_state._stepping;

        using geo_cxt_t =
            typename propagator_state_t::detector_type::geometry_context;
        const geo_cxt_t ctx{};

        inspector_state._nav_status.push_back(navigation.status());

        if (prop_state.param_type() == parameter_type::e_free) {
            return;
        }

        // Nothing has happened yet (first call of actor chain)
        if (stepping.path_length() < tol || stepping._s < tol) {
            return;
        }

        // Surface
        const auto sf = surface{*navigation.detector(),
                                stepping._bound_params.surface_link()};

        const auto free_vec =
            sf.bound_to_free_vector(ctx, stepping._bound_params.vector());

        const auto last_pos =
            detail::track_helper<matrix_operator>().pos(free_vec);

        free_track_parameters<algebra_t> free_params;
        free_params.set_vector(free_vec);

        const auto bvec =
            stepping._magnetic_field.at(last_pos[0], last_pos[1], last_pos[2]);
        const vector3 b{bvec[0], bvec[1], bvec[2]};

        detail::helix<algebra_t> hlx(free_params, &b);

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
GTEST_TEST(detray_propagator, propagator_line_stepper) {

    vecmem::host_memory_resource host_mr;
    toy_det_config toy_cfg{};
    toy_cfg.use_material_maps(false);
    const auto [d, names] = build_toy_detector(host_mr, toy_cfg);

    using navigator_t = navigator<decltype(d), navigation::print_inspector>;
    using stepper_t = line_stepper<algebra_t>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

    const point3 pos{0.f, 0.f, 0.f};
    const vector3 mom{1.f, 1.f, 0.f};
    free_track_parameters<algebra_t> track(pos, 0.f, mom, -1.f);

    propagator_t p{};

    propagator_t::state state(track, d);

    EXPECT_TRUE(p.propagate(state))
        << state._navigation.inspector().to_string() << std::endl;
}

/// Fixture for Runge-Kutta Propagation
class PropagatorWithRkStepper
    : public ::testing::TestWithParam<
          std::tuple<scalar_t, scalar_t, test::vector3>> {

    public:
    using generator_t =
        uniform_track_generator<free_track_parameters<algebra_t>>;

    /// Set the test environment up
    virtual void SetUp() {
        overstep_tol = std::get<0>(GetParam());
        step_constr = std::get<1>(GetParam());

        trk_gen_cfg.phi_steps(50u).theta_steps(50u);
        trk_gen_cfg.p_tot(10.f * unit<scalar_t>::GeV);
    }

    /// Clean up
    virtual void TearDown() {}

    protected:
    /// Detector configuration
    vecmem::host_memory_resource host_mr;

    /// Toy detector configuration
    toy_det_config toy_cfg = toy_det_config{}.n_brl_layers(4u).n_edc_layers(7u);

    /// Track generator configuration
    generator_t::configuration trk_gen_cfg{};

    /// Stepper configuration
    scalar_t overstep_tol, step_constr;
};

/// Test propagation in a constant magnetic field using a Runge-Kutta stepper
TEST_P(PropagatorWithRkStepper, rk4_propagator_const_bfield) {

    // Constant magnetic field type
    using bfield_t = bfield::const_field_t;

    // Toy detector
    using detector_t = detector<toy_metadata>;

    // Runge-Kutta propagation
    using navigator_t = navigator<detector_t /*, navigation::print_inspector*/>;
    using track_t = free_track_parameters<algebra_t>;
    using constraints_t = constrained_step<>;
    using policy_t = stepper_rk_policy;
    using stepper_t =
        rk_stepper<bfield_t::view_t, algebra_t, constraints_t, policy_t>;
    // Include helix actor to check track position/covariance
    using actor_chain_t =
        actor_chain<dtuple, helix_inspector, pathlimit_aborter,
                    parameter_transporter<algebra_t>,
                    pointwise_material_interactor<algebra_t>,
                    parameter_resetter<algebra_t>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Build detector
    toy_cfg.use_material_maps(false);
    toy_cfg.mapped_material(detray::vacuum<scalar_t>());
    const auto [det, names] = build_toy_detector(host_mr, toy_cfg);

    const bfield_t bfield = bfield::create_const_field(std::get<2>(GetParam()));

    // Propagator is built from the stepper and navigator
    propagation::config cfg{};
    cfg.navigation.overstep_tolerance = static_cast<float>(overstep_tol);
    propagator_t p{cfg};

    // Iterate through uniformly distributed momentum directions
    for (auto track : generator_t{trk_gen_cfg}) {

        assert(track.qop() != 0.f);

        // Generate second track state used for propagation with pathlimit
        track_t lim_track(track);

        // Build actor states: the helix inspector can be shared
        helix_inspector::state helix_insp_state{};
        pathlimit_aborter::state unlimted_aborter_state{};
        pathlimit_aborter::state pathlimit_aborter_state{path_limit};
        parameter_transporter<algebra_t>::state transporter_state{};
        pointwise_material_interactor<algebra_t>::state interactor_state{};
        parameter_resetter<algebra_t>::state resetter_state{};

        // Create actor states tuples
        auto actor_states =
            std::tie(helix_insp_state, unlimted_aborter_state,
                     transporter_state, interactor_state, resetter_state);
        auto sync_actor_states =
            std::tie(helix_insp_state, unlimted_aborter_state,
                     transporter_state, interactor_state, resetter_state);
        auto lim_actor_states =
            std::tie(helix_insp_state, pathlimit_aborter_state,
                     transporter_state, interactor_state, resetter_state);

        // Init propagator states
        propagator_t::state state(track, bfield, det);
        propagator_t::state lim_state(lim_track, bfield, det);

        // Set step constraints
        state._stepping.template set_constraint<step::constraint::e_accuracy>(
            step_constr);
        lim_state._stepping
            .template set_constraint<step::constraint::e_accuracy>(step_constr);

        // Propagate the entire detector
        state.do_debug = true;
        ASSERT_TRUE(p.propagate(state, actor_states))
            << state.debug_stream.str() << std::endl;
        //  << state._navigation.inspector().to_string() << std::endl;

        // Propagate with path limit
        ASSERT_NEAR(pathlimit_aborter_state.path_limit(), path_limit, tol);
        lim_state.do_debug = true;
        ASSERT_FALSE(p.propagate(lim_state, lim_actor_states))
            << lim_state.debug_stream.str() << std::endl;
        //<< lim_state._navigation.inspector().to_string() << std::endl;

        ASSERT_GE(std::abs(path_limit), lim_state._stepping._abs_path_length)
            << "Absolute path length: " << lim_state._stepping._abs_path_length
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

/// Test propagation in an inhomogenous magnetic field using a Runge-Kutta
/// stepper
TEST_P(PropagatorWithRkStepper, rk4_propagator_inhom_bfield) {

    // Magnetic field map using nearest neightbor interpolation
    using bfield_t = bfield::inhom_field_t;

    // Toy detector
    using detector_t = detector<toy_metadata>;

    // Runge-Kutta propagation
    using navigator_t = navigator<detector_t /*, navigation::print_inspector*/>;
    using track_t = free_track_parameters<algebra_t>;
    using constraints_t = constrained_step<>;
    using policy_t = stepper_rk_policy;
    using stepper_t =
        rk_stepper<bfield_t::view_t, algebra_t, constraints_t, policy_t>;
    // Include helix actor to check track position/covariance
    using actor_chain_t =
        actor_chain<dtuple, pathlimit_aborter, parameter_transporter<algebra_t>,
                    pointwise_material_interactor<algebra_t>,
                    parameter_resetter<algebra_t>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Build detector and magnetic field
    toy_cfg.use_material_maps(false);
    const auto [det, names] = build_toy_detector(host_mr, toy_cfg);
    const bfield_t bfield = bfield::create_inhom_field();

    // Propagator is built from the stepper and navigator
    propagation::config cfg{};
    cfg.navigation.overstep_tolerance = static_cast<float>(overstep_tol);
    propagator_t p{cfg};

    // Iterate through uniformly distributed momentum directions
    for (auto track : generator_t{trk_gen_cfg}) {
        // Genrate second track state used for propagation with pathlimit
        track_t lim_track(track);

        // Build actor states: the helix inspector can be shared
        pathlimit_aborter::state unlimted_aborter_state{};
        pathlimit_aborter::state pathlimit_aborter_state{path_limit};
        parameter_transporter<algebra_t>::state transporter_state{};
        pointwise_material_interactor<algebra_t>::state interactor_state{};
        parameter_resetter<algebra_t>::state resetter_state{};

        // Create actor states tuples
        auto actor_states = std::tie(unlimted_aborter_state, transporter_state,
                                     interactor_state, resetter_state);
        auto lim_actor_states =
            std::tie(pathlimit_aborter_state, transporter_state,
                     interactor_state, resetter_state);

        // Init propagator states
        propagator_t::state state(track, bfield, det);
        propagator_t::state lim_state(lim_track, bfield, det);

        // Set step constraints
        state._stepping.template set_constraint<step::constraint::e_accuracy>(
            step_constr);
        lim_state._stepping
            .template set_constraint<step::constraint::e_accuracy>(step_constr);

        // Propagate the entire detector
        state.do_debug = true;
        ASSERT_TRUE(p.propagate(state, actor_states))
            << state.debug_stream.str() << std::endl;
        //<< state._navigation.inspector().to_string() << std::endl;

        // Propagate with path limit
        ASSERT_NEAR(pathlimit_aborter_state.path_limit(), path_limit, tol);
        lim_state.do_debug = true;
        ASSERT_FALSE(p.propagate(lim_state, lim_actor_states))
            << lim_state.debug_stream.str() << std::endl;
        //<< lim_state._navigation.inspector().to_string() << std::endl;

        ASSERT_TRUE(lim_state._stepping.path_length() <
                    std::abs(path_limit) + tol)
            << "path length: " << lim_state._stepping.path_length()
            << ", path limit: " << path_limit << std::endl;
        //<< state._navigation.inspector().to_string() << std::endl;
    }
}

// No step size constraint
INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation1, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-0.05f * unit<scalar_t>::mm,
                                      std::numeric_limits<scalar_t>::max(),
                                      vector3{0.f * unit<scalar_t>::T,
                                              0.f * unit<scalar_t>::T,
                                              2.f * unit<scalar_t>::T})));

// Add some restrictions for more frequent navigation updates in the cases of
// non-z-aligned B-fields
INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation2, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-800.f * unit<scalar_t>::um,
                                      40.f * unit<scalar_t>::mm,
                                      vector3{0.f * unit<scalar_t>::T,
                                              1.f * unit<scalar_t>::T,
                                              1.f * unit<scalar_t>::T})));

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation3, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-800.f * unit<scalar_t>::um,
                                      40.f * unit<scalar_t>::mm,
                                      vector3{1.f * unit<scalar_t>::T,
                                              0.f * unit<scalar_t>::T,
                                              1.f * unit<scalar_t>::T})));

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation4, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-800.f * unit<scalar_t>::um,
                                      35.f * unit<scalar_t>::mm,
                                      vector3{1.f * unit<scalar_t>::T,
                                              1.f * unit<scalar_t>::T,
                                              1.f * unit<scalar_t>::T})));
