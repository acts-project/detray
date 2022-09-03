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

    /// Pathlimit for a single propagation workflow
    struct state {
        /// Absolute path limit
        scalar _path_limit = std::numeric_limits<scalar>::max();

        /// Set the path limit to a scalar @param pl
        DETRAY_HOST_DEVICE
        inline void set_path_limit(const scalar pl) { _path_limit = pl; }

        /// @returns this states remaining path length.
        DETRAY_HOST_DEVICE
        inline scalar path_limit() const { return _path_limit; }
    };

    /// Enforces the path limit on a stepper state
    ///
    /// @param abrt_state contains the path limit
    /// @param prop_state state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state &abrt_state,
                                       propagator_state_t &prop_state) const {
        auto &step_state = prop_state._stepping;
        auto &nav_state = prop_state._navigation;

        // Nothing left to do. Propagation will exit successfully
        if (nav_state.is_complete()) {
            return;
        }

        // Check the path limit
        abrt_state._path_limit -= std::abs(prop_state._stepping.step_size());
        if (abrt_state.path_limit() <= 0) {
            // printf("Abort: Stepper above maximal path length!\n");
            // Stop navigation
            prop_state._heartbeat &= nav_state.abort();
        }

        // Don't go over the path limit in the next step
        step_state.template set_constraint<step::constraint::e_aborter>(
            abrt_state.path_limit());
    }
};

struct looper : actor {

    struct state {
        scalar loop_length = 0;
        dindex target_surface_index = dindex_invalid;
    };

    /// Enforces the path limit on a stepper state
    ///
    /// @param abrt_state contains the path limit
    /// @param prop_state state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state &actor_state,
                                       propagator_state_t &prop_state) const {

        prop_state._heartbeat = true;
        auto &navigation = prop_state._navigation;
        auto &stepping = prop_state._stepping;
        navigation.set_full_trust();

        scalar residual = actor_state.loop_length - stepping.path_length();
        stepping.set_constraint(residual);

        typename propagator_state_t::navigator_state_type::intersection_t is;
        is.index = actor_state.target_surface_index;
        is.path = residual;
        auto &candidates = navigation.candidates();
        candidates.clear();
        candidates.push_back(is);
        navigation.set_next(candidates.begin());
        navigation.set_unknown();

        if (residual < std::abs(1e-5)) {
            prop_state._heartbeat = false;
            candidates.push_back({});
            navigation.next()++;
            navigation.set_on_module();
        }
    }
};

}  // namespace detray
