/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray definitions
#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"

// system includes
#include <climits>

namespace detray {

/// Struct that represents the most conservative navigation policy: alway re-
/// initialize the current volume
struct always_init : actor {

    struct state_type {};

    /// Sets the navigation trust level to 'no trust'
    ///
    /// @param pol_state not used
    /// @param propagation state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void operator()(
        const state_type & /*pol_state*/,
        propagator_state_t &propagation) const {
        propagation._navigation.set_no_trust();
    }
};

/// During guided navigation only the next surface should be re-evaluated. This
/// maps to the 'high trust' level in the navigator
struct guided_navigation : actor {

    struct state_type {};

    /// Sets the navigation trust level to 'no trust'
    ///
    /// @param pol_state not used
    /// @param propagation state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void operator()(
        const state_type & /*pol_state*/,
        propagator_state_t &propagation) const {
        propagation._navigation.set_high_trust();
    }
};

/// Default navigation update policy for the steppers: If a constraint has been
/// hit, lower the trustlevel to 'fair trust', otherwise stay in 'high trust'.
/// The reasoning is, that the track state might have changed much when a
/// constraint was triggered.
struct stepper_default_policy : actor {

    struct state_type {};

    /// Sets the navigation trust level depending on the step size limit
    ///
    /// @param pol_state not used
    /// @param propagation state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void operator()(
        state_type & /*pol_state*/, propagator_state_t &propagation) const {

        const auto &stepping = propagation._stepping;
        auto &navigation = propagation._navigation;

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

/// Specific navigation policy for the Runge-Kutta stepper: Use the relative
/// amount of step size correction as a measure for the change in direction of
/// the track state.
struct stepper_rk_correction_policy : actor {

    struct state_type {};

    /// Sets the navigation trust level depending on the step size limit
    ///
    /// @param pol_state not used
    /// @param propagation state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void operator()(
        state_type & /*pol_state*/, propagator_state_t &propagation) const {

        const auto &stepping = propagation._stepping;
        auto &navigation = propagation._navigation;

        scalar rel_correction =
            (stepping.step_size() - navigation()) / navigation();

        // Not a severe change to track state expected
        if (rel_correction > _threshold_no_trust) {
            // Re-evaluate only next candidate
            navigation.set_no_trust();
        }
        // Step size hit a constraint - the track state was probably changed a
        // lot
        else if (rel_correction > _threshold_fair_trust) {
            // Re-evaluate all candidates
            navigation.set_no_trust();
        } else {
            navigation.set_high_trust();
        }
    }

    scalar _threshold_fair_trust = 0.05;
    scalar _threshold_no_trust = 0.1;
};

}  // namespace detray