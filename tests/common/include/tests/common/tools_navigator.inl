/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/inspectors.hpp"

// GoogleTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <map>

namespace detray {

namespace {

// dummy propagator state
template <typename stepping_t, typename navigation_t>
struct prop_state {
    stepping_t _stepping;
    navigation_t _navigation;
};

/// Checks for a correct 'towards_surface' state
template <typename navigator_t, typename state_t = typename navigator_t::state>
inline void check_towards_surface(state_t &state, dindex vol_id,
                                  std::size_t n_candidates, dindex next_id) {
    ASSERT_EQ(state.status(), navigation::status::e_towards_object);
    ASSERT_EQ(state.volume(), vol_id);
    ASSERT_EQ(state.n_candidates(), n_candidates);
    // If we are towards some object, we have no current one (even if we are
    // geometrically still there)
    ASSERT_EQ(state.current_object(), dindex_invalid);
    // the portal is still the next object, since we did not step
    ASSERT_EQ(state.next_object(), next_id);
    ASSERT_TRUE((state.trust_level() == navigation::trust_level::e_full) or
                (state.trust_level() == navigation::trust_level::e_high));
}

/// Checks for a correct 'on_surface' state
template <typename navigator_t, typename state_t = typename navigator_t::state>
inline void check_on_surface(state_t &state, dindex vol_id,
                             std::size_t n_candidates, dindex current_id,
                             dindex next_id) {
    // The status is: on surface/towards surface if the next candidate is
    // immediately updated and set in the same update call
    ASSERT_TRUE(state.status() == navigation::status::e_on_module or
                state.status() == navigation::status::e_on_portal);
    // Points towards next candidate
    ASSERT_TRUE(std::abs(state()) > state.tolerance());
    ASSERT_EQ(state.volume(), vol_id);
    ASSERT_EQ(state.n_candidates(), n_candidates);
    ASSERT_EQ(state.current_object(), current_id);
    // points to the next surface now
    ASSERT_EQ(state.next_object(), next_id);
    ASSERT_EQ(state.trust_level(), navigation::trust_level::e_full);
}

/// Checks for a correctly handled volume switch
template <typename navigator_t, typename state_t = typename navigator_t::state>
inline void check_volume_switch(state_t &state, dindex vol_id) {
    // Switched to next volume
    ASSERT_EQ(state.volume(), vol_id);
    // The status is towards first surface in new volume
    ASSERT_EQ(state.status(), navigation::status::e_on_portal);
    // Kernel is newly initialized
    ASSERT_FALSE(state.is_exhausted());
    ASSERT_EQ(state.trust_level(), navigation::trust_level::e_full);
}

/// Checks an entire step onto the next surface
template <typename navigator_t, typename stepper_t, typename prop_state_t>
inline void check_step(navigator_t &nav, stepper_t &stepper,
                       prop_state_t &propagation, dindex vol_id,
                       std::size_t n_candidates, dindex current_id,
                       dindex next_id) {
    auto &navigation = propagation._navigation;

    // Step onto the surface in volume
    stepper.step(propagation);
    navigation.set_high_trust();
    // Stepper reduced trust level
    ASSERT_TRUE(navigation.trust_level() == navigation::trust_level::e_high);
    ASSERT_TRUE(nav.update(propagation));
    // Trust level is restored
    ASSERT_EQ(navigation.trust_level(), navigation::trust_level::e_full);
    // The status is on surface
    check_on_surface<navigator_t>(navigation, vol_id, n_candidates, current_id,
                                  next_id);
}

}  // anonymous namespace

}  // namespace detray

/// @note __plugin has to be defined with a preprocessor command

/// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, navigator) {
    using namespace detray;
    using namespace detray::navigation;
    using transform3 = __plugin::transform3<scalar>;

    vecmem::host_memory_resource host_mr;

    /// Tolerance for tests
    constexpr double tol = 0.01;

    std::size_t n_brl_layers = 4;
    std::size_t n_edc_layers = 3;
    auto toy_det = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);
    using detector_t = decltype(toy_det);
    using inspector_t = navigation::print_inspector;
    using navigator_t = navigator<detector_t, inspector_t>;
    using constraint_t = constrained_step<>;
    using stepper_t = line_stepper<transform3, constraint_t>;

    // test track
    point3 pos{0., 0., 0.};
    vector3 mom{1., 1., 0.};
    free_track_parameters<transform3> traj(pos, 0, mom, -1);

    stepper_t stepper;
    navigator_t nav(toy_det);

    prop_state<stepper_t::state, navigator_t::state> propagation{
        stepper_t::state{traj}, navigator_t::state{}};
    navigator_t::state &navigation = propagation._navigation;
    stepper_t::state &stepping = propagation._stepping;

    // Check that the state is unitialized
    // Default volume is zero
    ASSERT_EQ(navigation.volume(), 0u);
    // No surface candidates
    ASSERT_EQ(navigation.n_candidates(), 0u);
    // You can not trust the state
    ASSERT_EQ(navigation.trust_level(), trust_level::e_no_trust);
    // The status is unkown
    ASSERT_EQ(navigation.status(), status::e_unknown);

    //
    // beampipe
    //

    // Initialize navigation
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(nav.init(propagation));
    // The status is towards beampipe
    // Two candidates: beampipe and portal
    // First candidate is the beampipe
    check_towards_surface<navigator_t>(navigation, 0, 2, 0);
    // Distance to beampipe surface
    ASSERT_NEAR(navigation(), 19., tol);

    // Let's make half the step towards the beampipe
    stepping.template set_constraint<step::constraint::e_user>(navigation() *
                                                               0.5);
    stepper.step(propagation);
    // Navigation policy might reduce trust level to fair trust
    navigation.set_fair_trust();
    // Release user constraint again
    stepping.template release_step<step::constraint::e_user>();
    ASSERT_TRUE(navigation.trust_level() == trust_level::e_fair);
    // Re-navigate
    ASSERT_TRUE(nav.update(propagation));
    // Trust level is restored
    ASSERT_EQ(navigation.trust_level(), trust_level::e_full);
    // The status remains: towards surface
    check_towards_surface<navigator_t>(navigation, 0, 2, 0);
    // Distance to beampipe is now halved
    ASSERT_NEAR(navigation(), 9.5, tol);

    // Let's immediately update, nothing should change, as there is full trust
    ASSERT_TRUE(nav.update(propagation));
    check_towards_surface<navigator_t>(navigation, 0, 2, 0);
    ASSERT_NEAR(navigation(), 9.5, tol);

    // Now step onto the beampipe (idx 0)
    check_step(nav, stepper, propagation, 0, 2, 0, 7);
    // New target: Distance to the beampipe volume cylinder portal
    ASSERT_NEAR(navigation(), 8, tol);

    // Step onto portal 7 in volume 0
    stepper.step(propagation);
    navigation.set_high_trust();
    ASSERT_TRUE(navigation.trust_level() == trust_level::e_high);
    ASSERT_TRUE(nav.update(propagation));
    ASSERT_EQ(navigation.trust_level(), trust_level::e_full);

    //
    // barrel
    //

    // Last volume before we leave world
    dindex last_vol_id = 13;

    // maps volume id to the sequence of surfaces that the navigator encounters
    std::map<dindex, std::vector<dindex>> sf_sequences;

    // layer 1
    sf_sequences[7] = {594, 491, 475, 492, 476, 595};
    // gap 1
    sf_sequences[8] = {598, 599};
    // layer 2
    sf_sequences[9] = {1050, 845, 813, 846, 814, 1051};
    // gap 2
    sf_sequences[10] = {1054, 1055};
    // layer 3
    sf_sequences[11] = {1786, 1454, 1402, 1787};
    // gap 3
    sf_sequences[12] = {1790, 1791};
    // layer 4
    sf_sequences[last_vol_id] = {2886, 2388, 2310, 2887};

    // Every iteration steps through one barrel layer
    for (const auto &[vol_id, sf_seq] : sf_sequences) {
        // Includes the portal we are automatically on
        std::size_t n_candidates = sf_seq.size();

        // We switched to next barrel volume
        check_volume_switch<navigator_t>(navigation, vol_id);

        // The status is: on adjacent portal in volume, towards next candidate
        check_on_surface<navigator_t>(navigation, vol_id, n_candidates,
                                      sf_seq[0], sf_seq[1]);

        // Step through the module surfaces
        for (std::size_t sf = 1; sf < sf_seq.size() - 1; ++sf) {
            check_step(nav, stepper, propagation, vol_id, n_candidates,
                       sf_seq[sf], sf_seq[sf + 1]);
        }

        // Step onto the portal in volume
        stepper.step(propagation);
        navigation.set_high_trust();

        // Check agianst last volume
        if (vol_id == last_vol_id) {
            ASSERT_FALSE(nav.update(propagation));
            // The status is: exited
            ASSERT_EQ(navigation.status(), status::e_on_target);
            // Switch to next volume leads out of the detector world -> exit
            ASSERT_EQ(navigation.volume(), dindex_invalid);
            // We know we went out of the detector
            ASSERT_EQ(navigation.trust_level(), trust_level::e_full);
        } else {
            ASSERT_TRUE(nav.update(propagation));
        }
    }

    // Leave for debugging
    // std::cout << navigation.inspector().to_string() << std::endl;
    ASSERT_TRUE(navigation.is_complete()) << navigation.inspector().to_string();
}