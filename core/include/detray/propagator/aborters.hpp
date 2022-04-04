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
#include "detray/propagator/base_stepper.hpp"

namespace detray {

/// Aborter that checks whether the track has exceeded its pathlimit
struct pathlimit_aborter : actor {

    // Pathlimit for a single propagation workflow
    struct aborter_state {
        // Absolute path limit
        scalar _path_limit = std::numeric_limits<scalar>::max();

        /// Set the path limit to a scalar @param pl
        DETRAY_HOST_DEVICE
        inline void set_path_limit(const scalar pl) { _path_limit = pl; }

        /// @returns this states remaining path length.
        DETRAY_HOST_DEVICE
        inline scalar path_limit() const { return _path_limit; }
    };

    /// Broadcast state type to actor chain
    using state_type = aborter_state;

    /// Enforces the path limit on a stepper state
    ///
    /// @param abrt_state contains the path limit
    /// @param prop_state state of the propagation
    template <typename propagator_state_t>
    void operator()(const state_type &abrt_state,
                    propagator_state_t &prop_state) const {
        auto &step_state = prop_state._stepping;
        auto &nav_state = prop_state._navigation;

        // Nothing left to do. Propagation will exit successfully
        if (nav_state.is_complete()) {
            return;
        }
        // std::cout << "Path limit: " << abrt_state.path_limit() << ", track
        // length: "
        //<< std::abs(prop_state._stepping.path_length()) << std::endl;

        if (abrt_state.path_limit() <=
            std::abs(prop_state._stepping.path_length())) {
            printf("Abort: Stepper above maximal path length!\n");
            // Stop navigation
            prop_state._heartbeat &= nav_state.abort();
        }

        // Don't go over the path limit in the next step
        step_state.template set_constraint<step::constraint::e_aborter>(
            abrt_state.path_limit() -
            std::abs(prop_state._stepping.path_length()));
    }
};

}  // namespace detray