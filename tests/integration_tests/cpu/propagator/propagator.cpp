/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/propagator/propagator.hpp"

#include "detray/definitions/units.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/navigation/direct_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/tracks/trajectories.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/inspectors.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;
using vector3 = test::vector3;

namespace {

constexpr scalar tol{1e-3f};
constexpr scalar path_limit{5.f * unit<scalar>::cm};
constexpr std::size_t cache_size{navigation::default_cache_size};

/// Compare helical track positions for stepper
struct helix_inspector : actor {

    /// Keeps the state of a helix gun to calculate track positions
    struct state {
        /// navigation status for every step
        std::vector<navigation::status> _nav_status;
        scalar path_from_surface{0.f};
    };

    /// Check that the stepper remains on the right helical track for its pos.
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        state& inspector_state, const propagator_state_t& prop_state) const {

        const auto& navigation = prop_state._navigation;
        const auto& stepping = prop_state._stepping;

        typename propagator_state_t::detector_type::geometry_context ctx{};

        // Update inspector state
        inspector_state._nav_status.push_back(navigation.status());
        // The propagation does not start on a surface, skipp the inital path
        if (!stepping.bound_params().surface_link().is_invalid()) {
            inspector_state.path_from_surface += stepping.step_size();
        }

        // Nothing has happened yet (first call of actor chain)
        if (stepping.path_length() < tol ||
            inspector_state.path_from_surface < tol) {
            return;
        }

        if (stepping.bound_params().surface_link().is_invalid()) {
            return;
        }

        // Surface
        const auto sf = tracking_surface{
            navigation.detector(), stepping.bound_params().surface_link()};

        const free_track_parameters<test_algebra> free_params =
            sf.bound_to_free_vector(ctx, stepping.bound_params());

        const auto last_pos = free_params.pos();

        const auto bvec =
            stepping.field().at(last_pos[0], last_pos[1], last_pos[2]);
        const vector3 b{bvec[0], bvec[1], bvec[2]};

        detail::helix<test_algebra> hlx(free_params, b);

        const auto true_pos = hlx(inspector_state.path_from_surface);

        const point3 relative_error{1.f / inspector_state.path_from_surface *
                                    (stepping().pos() - true_pos)};

        ASSERT_NEAR(vector::norm(relative_error), 0.f, tol);

        auto true_J = hlx.jacobian(inspector_state.path_from_surface);

        for (unsigned int i = 0u; i < e_free_size; i++) {
            for (unsigned int j = 0u; j < e_free_size; j++) {
                ASSERT_NEAR(
                    getter::element(stepping.transport_jacobian(), i, j),
                    getter::element(true_J, i, j),
                    inspector_state.path_from_surface * tol * 10.f);
            }
        }
        // Reset path from surface
        if (navigation.is_on_sensitive() ||
            navigation.encountered_sf_material()) {
            inspector_state.path_from_surface = 0.f;
        }
    }
};

}  // anonymous namespace

/// Test basic functionality of the propagator using a straight line stepper
GTEST_TEST(detray_propagator, propagator_line_stepper) {

    vecmem::host_memory_resource host_mr;
    toy_det_config<scalar> toy_cfg{};
    toy_cfg.use_material_maps(false);
    const auto [d, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

    using navigator_t =
        caching_navigator<decltype(d), cache_size, navigation::print_inspector>;
    using stepper_t = line_stepper<test_algebra>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

    const point3 pos{0.f, 0.f, 0.f};
    const vector3 mom{1.f, 1.f, 0.f};
    free_track_parameters<test_algebra> track(pos, 0.f, mom, -1.f);

    propagation::config prop_cfg{};
    propagator_t p{prop_cfg};

    propagator_t::state state(track, d, prop_cfg.context);

    EXPECT_TRUE(p.propagate(state))
        << state._navigation.inspector().to_string() << std::endl;
}

/// Fixture for Runge-Kutta Propagation
class PropagatorWithRkStepper : public ::testing::TestWithParam<
                                    std::tuple<scalar, scalar, test::vector3>> {

    public:
    using generator_t =
        uniform_track_generator<free_track_parameters<test_algebra>>;

    /// Set the test environment up
    virtual void SetUp() {
        overstep_tol = std::get<0>(GetParam());
        step_constr = std::get<1>(GetParam());

        trk_gen_cfg.phi_steps(50u).theta_steps(50u);
        trk_gen_cfg.p_tot(10.f * unit<scalar>::GeV);
    }

    /// Clean up
    virtual void TearDown() { /* Do nothing */ }

    protected:
    /// Detector configuration
    vecmem::host_memory_resource host_mr;

    /// Toy detector configuration
    toy_det_config<scalar> toy_cfg =
        toy_det_config<scalar>{}.n_brl_layers(4u).n_edc_layers(7u);

    /// Track generator configuration
    generator_t::configuration trk_gen_cfg{};

    /// Stepper configuration
    scalar overstep_tol;
    scalar step_constr;
};

/// Test propagation in a constant magnetic field using a Runge-Kutta stepper
TEST_P(PropagatorWithRkStepper, rk4_propagator_const_bfield) {

    // Constant magnetic field type
    using bfield_t = bfield::const_field_t<scalar>;

    // Toy detector
    using detector_t = detector<test::toy_metadata>;

    // Runge-Kutta propagation
    using navigator_t =
        caching_navigator<detector_t, cache_size, navigation::print_inspector>;
    using track_t = free_track_parameters<test_algebra>;
    using constraints_t = constrained_step<scalar>;
    using policy_t = stepper_rk_policy<scalar>;
    using stepper_t =
        rk_stepper<bfield_t::view_t, test_algebra, constraints_t, policy_t>;
    // Include helix actor to check track position/covariance
    using actor_chain_t =
        actor_chain<helix_inspector, pathlimit_aborter<scalar>,
                    parameter_transporter<test_algebra>,
                    pointwise_material_interactor<test_algebra>,
                    parameter_resetter<test_algebra>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Build detector
    toy_cfg.use_material_maps(false);
    toy_cfg.mapped_material(detray::vacuum<scalar>());
    const auto [det, names] =
        build_toy_detector<test_algebra>(host_mr, toy_cfg);

    const bfield_t bfield = create_const_field<scalar>(std::get<2>(GetParam()));

    // Propagator is built from the stepper and navigator
    propagation::config cfg{};
    cfg.navigation.intersection.overstep_tolerance =
        static_cast<float>(overstep_tol);
    cfg.navigation.search_window = {3u, 3u};
    cfg.navigation.estimate_scattering_noise = false;
    propagator_t p{cfg};

    // Iterate through uniformly distributed momentum directions
    for (auto track : generator_t{trk_gen_cfg}) {

        assert(track.qop() != 0.f);

        // Generate second track state used for propagation with pathlimit
        track_t lim_track(track);

        // Build actor states: the helix inspector can be shared
        auto actor_states = actor_chain_t::make_default_actor_states();
        auto actor_states_lim = actor_chain_t::make_default_actor_states();
        auto actor_states_sync = actor_chain_t::make_default_actor_states();

        // Make sure the lim state is being terminated
        auto& pathlimit_aborter_state =
            detail::get<pathlimit_aborter<scalar>::state>(actor_states_lim);
        pathlimit_aborter_state.set_path_limit(path_limit);

        // Init propagator states
        propagator_t::state state(track, bfield, det);
        propagator_t::state sync_state(track, bfield, det);
        propagator_t::state lim_state(lim_track, bfield, det);

        state.do_debug = true;
        sync_state.do_debug = true;
        lim_state.do_debug = true;

        // Set step constraints
        state._stepping.template set_constraint<step::constraint::e_accuracy>(
            step_constr);
        sync_state._stepping
            .template set_constraint<step::constraint::e_accuracy>(step_constr);
        lim_state._stepping
            .template set_constraint<step::constraint::e_accuracy>(step_constr);

        // No multiple scattering simulated in this test
        using resetter_state_t = parameter_resetter<test_algebra>::state;
        detail::get<resetter_state_t>(actor_states).estimate_scattering_noise =
            false;
        detail::get<resetter_state_t>(actor_states_lim)
            .estimate_scattering_noise = false;

        // Propagate the entire detector
        ASSERT_TRUE(
            p.propagate(state, actor_chain_t::setup_actor_states(actor_states)))
            << state._navigation.inspector().to_string() << std::endl;

        // Propagate with path limit
        ASSERT_FALSE(p.propagate(
            lim_state, actor_chain_t::setup_actor_states(actor_states_lim)))
            << lim_state._navigation.inspector().to_string() << std::endl;

        ASSERT_GE(std::abs(path_limit), lim_state._stepping.abs_path_length())
            << "Absolute path length: " << lim_state._stepping.abs_path_length()
            << ", path limit: " << path_limit << std::endl;
        //<< state._navigation.inspector().to_string() << std::endl;
    }
}

/// Test propagation in an inhomogenous magnetic field using a Runge-Kutta
/// stepper
TEST_P(PropagatorWithRkStepper, rk4_propagator_inhom_bfield) {

    // Magnetic field map using nearest neightbor interpolation
    using bfield_t = bfield::inhom_field_t<scalar>;

    // Toy detector
    using detector_t = detector<test::toy_metadata>;

    // Runge-Kutta propagation
    using navigator_t =
        caching_navigator<detector_t, cache_size, navigation::print_inspector>;
    using track_t = free_track_parameters<test_algebra>;
    using constraints_t = constrained_step<scalar>;
    using policy_t = stepper_rk_policy<scalar>;
    using stepper_t =
        rk_stepper<bfield_t::view_t, test_algebra, constraints_t, policy_t>;
    // Include helix actor to check track position/covariance
    using actor_chain_t =
        actor_chain<pathlimit_aborter<scalar>,
                    parameter_transporter<test_algebra>,
                    pointwise_material_interactor<test_algebra>,
                    parameter_resetter<test_algebra>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Build detector and magnetic field
    toy_cfg.use_material_maps(false);
    const auto [det, names] =
        build_toy_detector<test_algebra>(host_mr, toy_cfg);
    const bfield_t bfield = create_inhom_field<scalar>();

    // Propagator is built from the stepper and navigator
    propagation::config cfg{};
    cfg.navigation.intersection.overstep_tolerance =
        static_cast<float>(overstep_tol);
    cfg.navigation.search_window = {3u, 3u};
    cfg.navigation.estimate_scattering_noise = false;
    propagator_t p{cfg};

    // Iterate through uniformly distributed momentum directions
    for (auto track : generator_t{trk_gen_cfg}) {
        // Genrate second track state used for propagation with pathlimit
        track_t lim_track(track);

        // Build actor states: the helix inspector can be shared
        pathlimit_aborter<scalar>::state unlimted_aborter_state{};
        pathlimit_aborter<scalar>::state pathlimit_aborter_state{path_limit};
        pointwise_material_interactor<test_algebra>::state interactor_state{};
        parameter_resetter<test_algebra>::state resetter_state{cfg};

        // Create actor states tuples
        auto actor_states = detray::tie(unlimted_aborter_state,
                                        interactor_state, resetter_state);
        auto lim_actor_states = detray::tie(pathlimit_aborter_state,
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
            << state._navigation.inspector().to_string() << std::endl;

        // Propagate with path limit
        ASSERT_NEAR(pathlimit_aborter_state.path_limit(), path_limit, tol);
        lim_state.do_debug = true;
        ASSERT_FALSE(p.propagate(lim_state, lim_actor_states))
            << lim_state._navigation.inspector().to_string() << std::endl;

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
    ::testing::Values(std::make_tuple(-100.f * unit<scalar>::um,
                                      std::numeric_limits<scalar>::max(),
                                      vector3{0.f * unit<scalar>::T,
                                              0.f * unit<scalar>::T,
                                              2.f * unit<scalar>::T})));

// Add some restrictions for more frequent navigation updates in the cases of
// non-z-aligned B-fields
INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation2, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-400.f * unit<scalar>::um,
                                      std::numeric_limits<scalar>::max(),
                                      vector3{0.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation3, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-400.f * unit<scalar>::um,
                                      std::numeric_limits<scalar>::max(),
                                      vector3{1.f * unit<scalar>::T,
                                              0.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_validation4, PropagatorWithRkStepper,
    ::testing::Values(std::make_tuple(-600.f * unit<scalar>::um,
                                      std::numeric_limits<scalar>::max(),
                                      vector3{1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T,
                                              1.f * unit<scalar>::T})));

/// Fixture for Runge-Kutta Propagation with direct navigator and toy detector
class PropagatorWithRkStepperDirectNavigatorToyDetector
    : public ::testing::TestWithParam<std::tuple<scalar, test::vector3>> {};

/// Test propagation in a constant magnetic field using a Runge-Kutta stepper
TEST_P(PropagatorWithRkStepperDirectNavigatorToyDetector, direct_navigator) {

    // Memory resource
    vecmem::host_memory_resource host_mr;

    // Track generator configuration
    using generator_t =
        uniform_track_generator<free_track_parameters<test_algebra>>;

    generator_t::configuration trk_gen_cfg{};
    trk_gen_cfg.phi_steps(50u).theta_steps(50u);
    trk_gen_cfg.p_tot(10.f * unit<scalar>::GeV);

    // Constant magnetic field type
    using bfield_t = bfield::const_field_t<scalar>;

    // Toy detector
    using detector_t = detector<test::toy_metadata>;
    using surface_t = typename detector_t::surface_type;

    // Runge-Kutta propagation
    using navigator_t =
        caching_navigator<detector_t, cache_size, navigation::print_inspector>;
    using constraints_t = constrained_step<scalar>;
    using policy_t = stepper_rk_policy<scalar>;
    using stepper_t =
        rk_stepper<bfield_t::view_t, test_algebra, constraints_t, policy_t>;
    // Include helix actor to check track position/covariance
    using actor_chain_t =
        actor_chain<parameter_transporter<test_algebra>,
                    pointwise_material_interactor<test_algebra>,
                    parameter_resetter<test_algebra>,
                    barcode_sequencer<surface_t>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    using direct_navigator_t = direct_navigator<detector_t>;
    using direct_propagator_t =
        propagator<stepper_t, direct_navigator_t, actor_chain_t>;

    // Build toy detector
    toy_det_config<scalar> toy_cfg =
        toy_det_config<scalar>{}.n_brl_layers(4u).n_edc_layers(7u);

    toy_cfg.use_material_maps(false);
    const auto [det, names] =
        build_toy_detector<test_algebra>(host_mr, toy_cfg);

    // Build mangetic field
    const bfield_t bfield = create_const_field<scalar>(std::get<1>(GetParam()));

    // Propagation configuration
    propagation::config cfg{};
    cfg.navigation.intersection.overstep_tolerance =
        static_cast<float>(std::get<0>(GetParam()));
    cfg.navigation.search_window = {3u, 3u};
    cfg.navigation.estimate_scattering_noise = false;
    propagation::config direct_cfg{};
    direct_cfg.navigation.intersection.min_mask_tolerance =
        std::numeric_limits<float>::max();
    direct_cfg.navigation.intersection.max_mask_tolerance =
        std::numeric_limits<float>::max();
    direct_cfg.navigation.intersection.overstep_tolerance =
        -std::numeric_limits<float>::max();
    propagator_t p{cfg};
    direct_propagator_t direct_p{direct_cfg};

    // Iterate through uniformly distributed momentum directions
    for (auto track : generator_t{trk_gen_cfg}) {

        // Build actor states: the helix inspector can be shared
        pointwise_material_interactor<test_algebra>::state interactor_state{};
        parameter_resetter<test_algebra>::state resetter_state{cfg};
        vecmem::data::vector_buffer<surface_t> seqs_buffer{
            100u, host_mr, vecmem::data::buffer_type::resizable};
        vecmem::data::vector_buffer<surface_t> seqs_forward_buffer{
            100u, host_mr, vecmem::data::buffer_type::resizable};
        vecmem::data::vector_buffer<surface_t> seqs_backward_buffer{
            100u, host_mr, vecmem::data::buffer_type::resizable};
        vecmem::copy m_copy;
        m_copy.setup(seqs_buffer)->wait();
        m_copy.setup(seqs_forward_buffer)->wait();
        m_copy.setup(seqs_backward_buffer)->wait();

        vecmem::device_vector<surface_t> seqs_device(seqs_buffer);
        vecmem::device_vector<surface_t> seqs_forward_device(
            seqs_forward_buffer);
        vecmem::device_vector<surface_t> seqs_backward_device(
            seqs_backward_buffer);

        barcode_sequencer<surface_t>::state sequencer_state(seqs_device);
        barcode_sequencer<surface_t>::state sequencer_forward_state(
            seqs_forward_device);
        barcode_sequencer<surface_t>::state sequencer_backward_state(
            seqs_backward_device);

        auto actor_states =
            detray::tie(interactor_state, sequencer_state, resetter_state);

        propagator_t::state state(track, bfield, det);
        navigator_t::state& navigation = state._navigation;

        // Propagate the entire detector
        // state.do_debug = true;
        ASSERT_TRUE(p.propagate(state, actor_states))
            << navigation.inspector().to_string();

        if (seqs_device.size() > 0) {

            auto direct_forward_actor_states = detray::tie(
                interactor_state, sequencer_forward_state, resetter_state);
            auto direct_backward_actor_states = detray::tie(
                interactor_state, sequencer_backward_state, resetter_state);

            direct_propagator_t::state direct_forward_state(track, bfield, det,
                                                            seqs_buffer);

            // direct_forward_state.do_debug = true;
            ASSERT_TRUE(direct_p.propagate(direct_forward_state,
                                           direct_forward_actor_states));

            // Check if all surfaces in the sequence are encountered
            ASSERT_TRUE(direct_forward_state._navigation.finished());
            ASSERT_EQ(sequencer_state.sequence().size(),
                      sequencer_forward_state.sequence().size());
            for (unsigned int i = 0; i < sequencer_state.sequence().size();
                 i++) {
                ASSERT_EQ(sequencer_state.sequence().at(i),
                          sequencer_forward_state.sequence().at(i));
            }

            const auto ptc = state._stepping.particle_hypothesis();
            ASSERT_EQ(
                ptc.pdg_num(),
                direct_forward_state._stepping.particle_hypothesis().pdg_num());
            const auto q = ptc.charge();

            // The initial momentum should be higher than the momentum at the
            // last surface
            ASSERT_TRUE(track.p(q) > state._stepping.bound_params().p(q));
            ASSERT_FLOAT_EQ(
                static_cast<float>(state._stepping.bound_params().p(q)),
                static_cast<float>(
                    direct_forward_state._stepping.bound_params().p(q)));

            direct_propagator_t::state direct_backward_state(
                direct_forward_state._stepping.bound_params(), bfield, det,
                seqs_buffer);
            direct_backward_state._navigation.set_direction(
                detray::navigation::direction::e_backward);

            ASSERT_TRUE(direct_p.propagate(direct_backward_state,
                                           direct_backward_actor_states));
            // Check if all surfaces in the sequence are encountered
            ASSERT_TRUE(direct_backward_state._navigation.finished());
            ASSERT_EQ(sequencer_state.sequence().size(),
                      sequencer_backward_state.sequence().size());
            for (unsigned int i = 0; i < sequencer_state.sequence().size();
                 i++) {
                unsigned int j = sequencer_state.sequence().size() - 1 - i;
                ASSERT_EQ(sequencer_state.sequence().at(i),
                          sequencer_backward_state.sequence().at(j));
            }

            ASSERT_NEAR(
                static_cast<float>(track.p(q)),
                static_cast<float>(
                    direct_backward_state._stepping.bound_params().p(q)),
                static_cast<float>(track.p(q)) * 0.0002f);
            ASSERT_TRUE(direct_backward_state._stepping.bound_params().p(q) >
                        direct_forward_state._stepping.bound_params().p(q));
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_direct_navigator_toy_detector,
    PropagatorWithRkStepperDirectNavigatorToyDetector,
    ::testing::Values(
        std::make_tuple(-100.f * unit<scalar>::um,
                        vector3{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                2.f * unit<scalar>::T}),
        std::make_tuple(-400.f * unit<scalar>::um,
                        vector3{0.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                1.f * unit<scalar>::T}),
        std::make_tuple(-400.f * unit<scalar>::um,
                        vector3{1.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                1.f * unit<scalar>::T}),
        std::make_tuple(-600.f * unit<scalar>::um,
                        vector3{1.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                1.f * unit<scalar>::T})));

/// Fixture for Runge-Kutta Propagation with direct navigator and wire chamber
class PropagatorWithRkStepperDirectNavigatorWireChamber
    : public ::testing::TestWithParam<std::tuple<scalar, test::vector3>> {};

/// Test propagation in a constant magnetic field using a Runge-Kutta stepper
TEST_P(PropagatorWithRkStepperDirectNavigatorWireChamber, direct_navigator) {

    // Memory resource
    vecmem::host_memory_resource host_mr;

    // Track generator configuration
    using generator_t =
        uniform_track_generator<free_track_parameters<test_algebra>>;

    generator_t::configuration trk_gen_cfg{};
    trk_gen_cfg.phi_steps(50u).theta_steps(50u);
    trk_gen_cfg.p_tot(10.f * unit<scalar>::GeV);

    // Constant magnetic field type
    using bfield_t = bfield::const_field_t<scalar>;

    // Toy detector
    using detector_t = detector<test::default_metadata>;
    using surface_t = typename detector_t::surface_type;

    // Runge-Kutta propagation
    using navigator_t =
        caching_navigator<detector_t, cache_size, navigation::print_inspector>;
    using constraints_t = constrained_step<scalar>;
    using policy_t = stepper_rk_policy<scalar>;
    using stepper_t =
        rk_stepper<bfield_t::view_t, test_algebra, constraints_t, policy_t>;
    // Include helix actor to check track position/covariance
    using actor_chain_t =
        actor_chain<parameter_transporter<test_algebra>,
                    pointwise_material_interactor<test_algebra>,
                    parameter_resetter<test_algebra>,
                    barcode_sequencer<surface_t>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    using direct_navigator_t = direct_navigator<detector_t>;
    using direct_propagator_t =
        propagator<stepper_t, direct_navigator_t, actor_chain_t>;

    // Build wire chamber
    wire_chamber_config<scalar> wire_cfg{};
    auto [det, names] = build_wire_chamber<test::algebra>(host_mr, wire_cfg);

    // Build mangetic field
    const bfield_t bfield = create_const_field<scalar>(std::get<1>(GetParam()));

    // Propagation configuration
    propagation::config cfg{};
    cfg.navigation.intersection.overstep_tolerance =
        static_cast<float>(std::get<0>(GetParam()));
    cfg.navigation.search_window = {5u, 5u};
    cfg.navigation.estimate_scattering_noise = false;
    propagation::config direct_cfg{};
    direct_cfg.navigation.intersection.min_mask_tolerance =
        std::numeric_limits<float>::max();
    direct_cfg.navigation.intersection.max_mask_tolerance =
        std::numeric_limits<float>::max();
    direct_cfg.navigation.intersection.overstep_tolerance =
        -std::numeric_limits<float>::max();
    propagator_t p{cfg};
    direct_propagator_t direct_p{direct_cfg};

    // Iterate through uniformly distributed momentum directions
    for (auto track : generator_t{trk_gen_cfg}) {

        // Build actor states: the helix inspector can be shared
        pointwise_material_interactor<test_algebra>::state interactor_state{};
        parameter_resetter<test_algebra>::state resetter_state{cfg};
        vecmem::data::vector_buffer<surface_t> seqs_buffer{
            100u, host_mr, vecmem::data::buffer_type::resizable};
        vecmem::data::vector_buffer<surface_t> seqs_forward_buffer{
            100u, host_mr, vecmem::data::buffer_type::resizable};
        vecmem::data::vector_buffer<surface_t> seqs_backward_buffer{
            100u, host_mr, vecmem::data::buffer_type::resizable};
        vecmem::copy m_copy;
        m_copy.setup(seqs_buffer)->wait();
        m_copy.setup(seqs_forward_buffer)->wait();
        m_copy.setup(seqs_backward_buffer)->wait();

        vecmem::device_vector<surface_t> seqs_device(seqs_buffer);
        vecmem::device_vector<surface_t> seqs_forward_device(
            seqs_forward_buffer);
        vecmem::device_vector<surface_t> seqs_backward_device(
            seqs_backward_buffer);

        barcode_sequencer<surface_t>::state sequencer_state(seqs_device);
        barcode_sequencer<surface_t>::state sequencer_forward_state(
            seqs_forward_device);
        barcode_sequencer<surface_t>::state sequencer_backward_state(
            seqs_backward_device);

        auto actor_states =
            detray::tie(interactor_state, sequencer_state, resetter_state);

        propagator_t::state state(track, bfield, det);
        navigator_t::state& navigation = state._navigation;

        // Propagate the entire detector
        ASSERT_TRUE(p.propagate(state, actor_states))
            << navigation.inspector().to_string();

        if (seqs_device.size() > 0) {

            auto direct_forward_actor_states = detray::tie(
                interactor_state, sequencer_forward_state, resetter_state);
            auto direct_backward_actor_states = detray::tie(
                interactor_state, sequencer_backward_state, resetter_state);

            direct_propagator_t::state direct_forward_state(track, bfield, det,
                                                            seqs_buffer);

            // direct_forward_state.do_debug = true;
            ASSERT_TRUE(direct_p.propagate(direct_forward_state,
                                           direct_forward_actor_states));

            // Check if all surfaces in the sequence are encountered
            ASSERT_TRUE(direct_forward_state._navigation.finished());
            ASSERT_EQ(sequencer_state.sequence().size(),
                      sequencer_forward_state.sequence().size());
            for (unsigned int i = 0; i < sequencer_state.sequence().size();
                 i++) {
                ASSERT_EQ(sequencer_state.sequence().at(i),
                          sequencer_forward_state.sequence().at(i));
            }

            const auto ptc = state._stepping.particle_hypothesis();
            ASSERT_EQ(
                ptc.pdg_num(),
                direct_forward_state._stepping.particle_hypothesis().pdg_num());
            const auto q = ptc.charge();

            // The initial momentum should be higher than or equal to the
            // momentum at the last surface
            ASSERT_GE(track.p(q), state._stepping.bound_params().p(q));
            ASSERT_NEAR(
                static_cast<float>(state._stepping.bound_params().p(q)),
                static_cast<float>(
                    direct_forward_state._stepping.bound_params().p(q)),
                static_cast<float>(state._stepping.bound_params().p(q)) *
                    2e-5f);

            direct_propagator_t::state direct_backward_state(
                direct_forward_state._stepping.bound_params(), bfield, det,
                seqs_buffer);
            direct_backward_state._navigation.set_direction(
                detray::navigation::direction::e_backward);

            ASSERT_TRUE(direct_p.propagate(direct_backward_state,
                                           direct_backward_actor_states));

            // Check if all surfaces in the sequence are encountered
            ASSERT_TRUE(direct_backward_state._navigation.finished());
            ASSERT_EQ(sequencer_state.sequence().size(),
                      sequencer_backward_state.sequence().size());
            for (unsigned int i = 0; i < sequencer_state.sequence().size();
                 i++) {
                unsigned int j = sequencer_state.sequence().size() - 1 - i;
                EXPECT_EQ(sequencer_state.sequence().at(i),
                          sequencer_backward_state.sequence().at(j));
            }

            ASSERT_NEAR(
                static_cast<float>(track.p(q)),
                static_cast<float>(
                    direct_backward_state._stepping.bound_params().p(q)),
                static_cast<float>(track.p(q)) * 0.0002f);
            ASSERT_GE(direct_backward_state._stepping.bound_params().p(q),
                      direct_forward_state._stepping.bound_params().p(q));
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_direct_navigator_wire_chamber,
    PropagatorWithRkStepperDirectNavigatorWireChamber,
    ::testing::Values(
        std::make_tuple(-100.f * unit<scalar>::um,
                        vector3{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                2.f * unit<scalar>::T}),
        std::make_tuple(-400.f * unit<scalar>::um,
                        vector3{0.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                1.f * unit<scalar>::T}),
        std::make_tuple(-400.f * unit<scalar>::um,
                        vector3{1.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                1.f * unit<scalar>::T}),
        std::make_tuple(-600.f * unit<scalar>::um,
                        vector3{1.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                1.f * unit<scalar>::T})));
