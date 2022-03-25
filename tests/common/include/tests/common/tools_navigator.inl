/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <map>

#include "detray/core/mask_store.hpp"
#include "detray/tools/line_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "detray/tools/track.hpp"
#include "detray/utils/indexing.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/inspectors.hpp"

namespace detray {

/** Checks for a correct 'towards_surface' state */
template <typename navigator_t, typename state_t = typename navigator_t::state>
inline void check_towards_surface(state_t &state, dindex vol_id,
                                  std::size_t n_candidates, dindex next_id) {
    ASSERT_EQ(state.status(), navigation::status::e_towards_object);
    ASSERT_EQ(state.volume(), vol_id);
    ASSERT_EQ(state.candidates().size(), n_candidates);
    // If we are towards some object, we have no current one (even if we are
    // geometrically still there)
    ASSERT_EQ(state.current_object(), dindex_invalid);
    // the portal is still the next object, since we did not step
    ASSERT_EQ(state.next()->index, next_id);
    ASSERT_TRUE((state.trust_level() == navigation::trust_level::e_full) or
                (state.trust_level() == navigation::trust_level::e_high));
}

/** Checks for a correct 'on_surface' state */
template <typename navigator_t, typename state_t = typename navigator_t::state>
inline void check_on_surface(state_t &state, dindex vol_id,
                             std::size_t n_candidates, dindex /*current_id*/,
                             dindex next_id) {
    // The status is: on surface/towards surface if the next candidate is
    // immediately updated and set in the same update call
    ASSERT_EQ(state.status(), navigation::status::e_towards_object);
    // Points towards next candidate
    ASSERT_TRUE(std::abs(state()) > state.tolerance());
    ASSERT_EQ(state.volume(), vol_id);
    ASSERT_EQ(state.candidates().size(), n_candidates);
    ASSERT_EQ(state.current_object(), dindex_invalid /*current_id*/);
    // points to the next surface now
    ASSERT_EQ(state.next()->index, next_id);
    ASSERT_EQ(state.trust_level(), navigation::trust_level::e_full);
}

/** Checks for a correctly handled volume switch */
template <typename navigator_t, typename state_t = typename navigator_t::state>
inline void check_volume_switch(state_t &state, dindex vol_id) {
    // Switched to next volume
    ASSERT_EQ(state.volume(), vol_id);
    // The status is towards first surface in new volume
    ASSERT_EQ(state.status(), navigation::status::e_towards_object);
    // Kernel is newly initialized
    ASSERT_FALSE(state.is_exhausted());
    ASSERT_EQ(state.trust_level(), navigation::trust_level::e_full);
}

/** Checks an entire step onto the next surface */
template <typename navigator_t, typename stepper_t, typename nav_state_t,
          typename stepper_state_t>
inline void check_step(navigator_t &n, stepper_t &s, nav_state_t &n_state,
                       stepper_state_t &s_state, dindex vol_id,
                       std::size_t n_candidates, dindex current_id,
                       dindex next_id) {

    // Step onto the surface in volume
    s.step(s_state, n_state);
    // Stepper reduced trust level
    ASSERT_TRUE(n_state.trust_level() == navigation::trust_level::e_high);
    ASSERT_TRUE(n.update(n_state, s_state));
    // Trust level is restored
    ASSERT_EQ(n_state.trust_level(), navigation::trust_level::e_full);
    // The status is on surface
    check_on_surface<navigator_t>(n_state, vol_id, n_candidates, current_id,
                                  next_id);
}

}  // namespace detray

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, navigator) {
    using namespace detray;
    using namespace detray::navigation;

    vecmem::host_memory_resource host_mr;

    /** Tolerance for tests */
    constexpr double tol = 0.01;

    std::size_t n_brl_layers = 4;
    std::size_t n_edc_layers = 3;
    auto toy_det = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);
    using detector_t = decltype(toy_det);
    // using inspector_t = navigation::void_inspector;
    using inspector_t = navigation::print_inspector;
    using navigator_t = navigator<detector_t, inspector_t>;
    navigator_t n(toy_det);
    using stepper_t = line_stepper<free_track_parameters>;

    // test track
    point3 pos{0., 0., 0.};
    vector3 mom{1., 1., 0.};
    free_track_parameters traj(pos, 0, mom, -1);
    stepper_t s;
    typename stepper_t::state s_state(traj);

    navigator_t::state n_state;

    // Check that the state is unitialized
    // Default volume is zero
    ASSERT_EQ(n_state.volume(), 0u);
    // No surface candidates
    ASSERT_EQ(n_state.candidates().size(), 0u);
    // You can not trust the state
    ASSERT_EQ(n_state.trust_level(), trust_level::e_no_trust);
    // The status is unkown
    ASSERT_EQ(n_state.status(), status::e_unknown);

    //
    // beampipe
    //

    // Initialize navigation
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(n.init(n_state, s_state));
    // The status is towards beampipe
    // Two candidates: beampipe and portal
    // First candidate is the beampipe
    check_towards_surface<navigator_t>(n_state, 0, 2, 0);
    // Distance to beampipe surface
    ASSERT_NEAR(n_state(), 19., tol);

    // Let's immediately target, nothing should change, as there is full trust
    ASSERT_TRUE(n.update(n_state, s_state));
    check_towards_surface<navigator_t>(n_state, 0, 2, 0);
    ASSERT_NEAR(n_state(), 19., tol);

    // Let's make half the step towards the beampipe
    s.step(s_state, n_state, n_state() * 0.5);
    // Stepper reduced trust level (hit step constrint -> only fair trust)
    ASSERT_TRUE(n_state.trust_level() == trust_level::e_fair);
    // Re-navigate
    ASSERT_TRUE(n.update(n_state, s_state));
    // Trust level is restored
    ASSERT_EQ(n_state.trust_level(), trust_level::e_full);
    // The status remains: towards surface
    check_towards_surface<navigator_t>(n_state, 0, 2, 0);
    // Distance to beampipe is now halved
    ASSERT_NEAR(n_state(), 9.5, tol);

    // Let's immediately target, nothing should change, as there is full trust
    ASSERT_TRUE(n.update(n_state, s_state));
    check_towards_surface<navigator_t>(n_state, 0, 2, 0);
    ASSERT_NEAR(n_state(), 9.5, tol);

    // Now step onto the beampipe (idx 0)
    check_step(n, s, n_state, s_state, 0, 2, 0, 7);
    // New target: Distance to the beampipe volume cylinder portal
    ASSERT_NEAR(n_state(), 8, tol);

    // Step onto portal 7 in volume 0
    s.step(s_state, n_state);
    ASSERT_TRUE(n_state.trust_level() == trust_level::e_high);
    ASSERT_TRUE(n.update(n_state, s_state));
    ASSERT_EQ(n_state.trust_level(), trust_level::e_full);

    //
    // barrel
    //

    // Last volume before we leave world
    dindex last_vol_id = 13;

    // maps volume id to the sequence of surfaces that the navigator encounters
    std::map<dindex, std::vector<dindex>> sf_sequences;

    // layer 1
    sf_sequences[7] = {491, 475, 492, 476, 595};
    // gap 1
    sf_sequences[8] = {599};
    // layer 2
    sf_sequences[9] = {845, 813, 846, 814, 1051};
    // gap 2
    sf_sequences[10] = {1055};
    // layer 3
    sf_sequences[11] = {1454, 1402, 1787};
    // gap 3
    sf_sequences[12] = {1791};
    // layer 4
    sf_sequences[last_vol_id] = {2388, 2310, 2887};

    // Every iteration steps through one barrel layer
    for (const auto &[vol_id, sf_seq] : sf_sequences) {
        // Includes the portal we are automatically on
        std::size_t n_candidates = sf_seq.size() + 1;

        // We switched to next barrel volume
        check_volume_switch<navigator_t>(n_state, vol_id);

        // The status is: on adjacent portal in volume, towards next candidate
        // This includes overlapping modules and the adjacent portal we are
        // already on
        check_towards_surface<navigator_t>(n_state, vol_id, n_candidates,
                                           sf_seq.front());

        // Step through the module surfaces
        for (std::size_t sf = 0; sf < sf_seq.size() - 1; ++sf) {
            check_step(n, s, n_state, s_state, vol_id, n_candidates, sf_seq[sf],
                       sf_seq[sf + 1]);
        }

        // Step onto the portal in volume
        s.step(s_state, n_state);

        // Check agianst last volume
        if (vol_id == last_vol_id) {
            ASSERT_FALSE(n.update(n_state, s_state));
            // The status is: exited
            ASSERT_EQ(n_state.status(), status::e_exit);
            // Switch to next volume leads out of the detector world -> exit
            ASSERT_EQ(n_state.volume(), dindex_invalid);
            // We know we went out of the detector
            ASSERT_EQ(n_state.trust_level(), trust_level::e_full);
        } else {
            ASSERT_TRUE(n.update(n_state, s_state));
        }
    }

    // std::cout << n_state.inspector().to_string() << std::endl;
}
