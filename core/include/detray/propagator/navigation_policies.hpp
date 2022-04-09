/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray definitions
#include <climits>

#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"

namespace detray {

/// Struct that represents the most conservative navigation policy: alway re-
/// initialize the current volume
struct always_init : actor {

    struct state_type {};

    /// Sets the navigation trust level to 'no trust'
    ///
    /// @param pol_state not used
    /// @param prop_state state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(const state_type & /*pol_state*/,
                                       propagator_state_t &prop_state) const {
        prop_state._navigation.set_no_trust();
    }
};

/// During guided navigation only the next surface should be re-evaluated. This
/// maps to the 'high trust' level in the navigator
struct guided_navigation : actor {

    struct state_type {};

    /// Sets the navigation trust level to 'no trust'
    ///
    /// @param pol_state not used
    /// @param prop_state state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(const state_type & /*pol_state*/,
                                       propagator_state_t &prop_state) const {
        prop_state._navigation.set_high_trust();
    }
};

namespace step {

/// Default navigation update policy for the steppers: If a constraint has been
/// hit, lower the trustlevel to 'fair trust', otherwise stay in 'high trust'.
/// The reasoning is, that the track state might have changed much when a
/// constraint was triggered.
struct default_policy : actor {

    /// State of the policy keeps track of the path length
    struct policy_state {
        scalar _path_length{0};
        scalar _tolerance{std::numeric_limits<scalar>::epsilon()};
    };

    using state_type = policy_state;

    /// Sets the navigation trust level to 'no trust'
    ///
    /// @param pol_state not used
    /// @param prop_state state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void operator()(
        state_type &pol_state, propagator_state_t &prop_state) const {

        const auto &stepping = prop_state._stepping;

        // This is run after the navigator update - do nothing
        if (std::abs(pol_state._path_length - stepping.path_length()) <
            pol_state._tolerance) {
            return;
        }

        // Update path length when the stepper advanced the track state
        pol_state._path_length = stepping.path_length();

        auto &navigation = prop_state._navigation;

        // Not a severe change to track state expected
        if (std::abs(stepping.step_size()) <
            std::abs(
                stepping.constraints().template size<>(stepping.direction()))) {
            // Re-evaluate only next candidate
            navigation.set_high_trust();
        }
        // Step size hit a constraint - the track state was probably changed a
        // lot
        else {
            // Re-evaluate all candidates
            navigation.set_fair_trust();
        }
    }
};

}  // namespace step

}  // namespace detray