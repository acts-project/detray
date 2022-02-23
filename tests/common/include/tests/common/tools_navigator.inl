/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <map>

#include "detray/core/mask_store.hpp"
#include "detray/tools/base_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "detray/tools/track.hpp"
#include "detray/utils/indexing.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

namespace detray {

/** Checks for a correct 'towards_surface' state */
template <typename navigator_t, typename state_t = typename navigator_t::state>
inline void check_towards_surface(state_t &state, dindex vol_id,
                                  std::size_t n_candidates, dindex next_id) {
    ASSERT_EQ(state.status(), navigator_t::navigation_status::e_towards_object);
    ASSERT_EQ(state.volume(), vol_id);
    ASSERT_EQ(state.candidates().size(), n_candidates);
    // If we are towards some object, we have no current one (even if we are
    // geometrically still there)
    ASSERT_EQ(state.on_object(), dindex_invalid);
    // the portal is still the next object, since we did not step
    ASSERT_EQ(state.next()->index, next_id);
    ASSERT_TRUE((state.trust_level() ==
                 navigator_t::navigation_trust_level::e_full_trust) or
                (state.trust_level() ==
                 navigator_t::navigation_trust_level::e_high_trust));
}

/** Checks for a correct 'on_surface' state */
template <typename navigator_t, typename state_t = typename navigator_t::state>
inline void check_on_surface(state_t &state, dindex vol_id,
                             std::size_t n_candidates, dindex current_id,
                             dindex next_id) {
    // The status is: on surface
    ASSERT_EQ(state.status(), navigator_t::navigation_status::e_on_object);
    ASSERT_TRUE(std::abs(state()) < state.tolerance());
    ASSERT_EQ(state.volume(), vol_id);
    ASSERT_EQ(state.candidates().size(), n_candidates);
    ASSERT_EQ(state.on_object(), current_id);
    // points to the next surface now
    ASSERT_EQ(state.next()->index, next_id);
    ASSERT_EQ(state.trust_level(),
              navigator_t::navigation_trust_level::e_high_trust);
}

/** Checks for a correctly handled volume switch */
template <typename navigator_t, typename state_t = typename navigator_t::state>
inline void check_volume_switch(state_t &state, dindex vol_id) {
    // The status is on portal
    ASSERT_EQ(state.status(), navigator_t::navigation_status::e_on_object);
    // Switched to next volume
    ASSERT_EQ(state.volume(), vol_id);
    // Kernel is exhaused, and trust level is gone
    ASSERT_EQ(state.next(), state.candidates().end());
    ASSERT_EQ(state.trust_level(),
              navigator_t::navigation_trust_level::e_no_trust);
}

/** Checks an entire step to next barrel surface */
template <typename navigator_t, typename state_t, typename stepper_state_t>
inline void check_step(navigator_t &nav, state_t &state,
                       stepper_state_t &stepping, dindex vol_id,
                       std::size_t n_candidates, dindex current_id,
                       dindex next_id) {

    auto &trck = stepping();

    // Step onto the surface in volume
    trck.set_pos(trck.pos() + state() * trck.dir());
    state.set_trust_level(navigator_t::navigation_trust_level::e_high_trust);
    ASSERT_TRUE(nav.status(state, stepping));
    // The status is: on surface 491
    check_on_surface<navigator_t>(state, vol_id, n_candidates, current_id,
                                  next_id);
    ASSERT_EQ(state.trust_level(),
              navigator_t::navigation_trust_level::e_high_trust);

    // Let's target - i.e. update the distance to next_id
    ASSERT_TRUE(nav.target(state, stepping));
    // Should be on our way to the next ovelapping module
    // this is still the next surface, since we did not step
    check_towards_surface<navigator_t>(state, vol_id, n_candidates, next_id);
}

}  // namespace detray

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, single_type_navigator) {
    using namespace detray;
    vecmem::host_memory_resource host_mr;

    /** Tolerance for tests */
    constexpr double tol = 0.01;

    std::size_t n_brl_layers = 4;
    std::size_t n_edc_layers = 3;
    auto toy_det = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);
    navigator n(toy_det);
    using toy_navigator = decltype(n);
    using stepper = base_stepper<free_track_parameters>;

    // test track
    point3 pos{0., 0., 0.};
    vector3 mom{1., 1., 0.};
    free_track_parameters traj(pos, 0, mom, -1);
    typename stepper::state stepping(traj);

    toy_navigator::state state;
    // Set initial volume (no grid yet)
    state.set_volume(0u);

    // Check that the state is unitialized
    // Volume is invalid
    ASSERT_EQ(state.volume(), 0u);
    // No surface candidates
    ASSERT_EQ(state.candidates().size(), 0u);
    // You can not trust the state
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_no_trust);
    // The status is unkown
    ASSERT_EQ(state.status(), toy_navigator::navigation_status::e_unknown);

    //
    // beampipe
    //

    // Initial status call
    bool heartbeat = n.status(state, stepping);
    // Test that the navigator has a heartbeat
    ASSERT_TRUE(heartbeat);
    // The status is towards beampipe
    // Two candidates: beampipe and portal
    // First candidate is the beampipe
    check_towards_surface<toy_navigator>(state, 0, 2, 0);
    ASSERT_NEAR(state(), 19., tol);

    // Let's immediately target, nothing should change, as there is full trust
    heartbeat = n.target(state, stepping);
    ASSERT_TRUE(heartbeat);
    // The status remains: towards surface
    check_towards_surface<toy_navigator>(state, 0, 2, 0);
    ASSERT_NEAR(state(), 19., tol);

    // Let's make half the step towards the beampipe
    traj.set_pos(traj.pos() + 0.5 * state() * traj.dir());
    // Could be externally set by actor (in the future)
    state.set_trust_level(toy_navigator::navigation_trust_level::e_high_trust);
    ASSERT_TRUE(n.status(state, stepping));
    // The status remains: towards surface
    check_towards_surface<toy_navigator>(state, 0, 2, 0);
    ASSERT_NEAR(state(), 9.5, tol);
    // Trust level is restored
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Let's immediately target, nothing should change, as there is full trust
    ASSERT_TRUE(n.target(state, stepping));
    check_towards_surface<toy_navigator>(state, 0, 2, 0);
    ASSERT_NEAR(state(), 9.5, tol);
    // Trust level is restored
    ASSERT_EQ(state.trust_level(),
              toy_navigator::navigation_trust_level::e_full_trust);

    // Now step onto the beampipe (idx 0)
    check_step(n, state, stepping, 0, 2, 0, 7);

    // Step onto portal 7 in volume 0
    traj.set_pos(traj.pos() + state() * traj.dir());
    ASSERT_TRUE(n.status(state, stepping));
    // We are in the first barrel layer now (vol 7)
    check_volume_switch<toy_navigator>(state, 7);

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
        check_volume_switch<toy_navigator>(state, vol_id);

        // New target call will initialize volume
        ASSERT_TRUE(n.target(state, stepping));
        // The status is: on adjacent portal in volume, towards next candidate
        // This includes overlapping modules and the adjacent portal we are
        // already on
        check_towards_surface<toy_navigator>(state, vol_id, n_candidates,
                                             sf_seq.front());

        // Step through the module surfaces
        for (std::size_t sf = 0; sf < sf_seq.size() - 1; ++sf) {
            check_step(n, state, stepping, vol_id, n_candidates, sf_seq[sf],
                       sf_seq[sf + 1]);
        }

        // Step onto the portal in volume
        traj.set_pos(traj.pos() + state() * traj.dir());
        state.set_trust_level(
            toy_navigator::navigation_trust_level::e_high_trust);

        // Check agianst last volume
        if (vol_id != last_vol_id) {
            ASSERT_TRUE(n.status(state, stepping));
        }
        // we leave the detector
        else {
            ASSERT_FALSE(n.status(state, stepping));
            // The status is: on portal
            ASSERT_EQ(state.status(),
                      toy_navigator::navigation_status::e_on_target);
            // Switch to next volume leads out of the detector world -> exit
            ASSERT_EQ(state.volume(), dindex_invalid);
            // We know we went out of the detector
            ASSERT_EQ(state.trust_level(),
                      toy_navigator::navigation_trust_level::e_full_trust);
        }
    }
}
