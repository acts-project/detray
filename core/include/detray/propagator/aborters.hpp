/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray definitions
#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/actor_chain.hpp"

#include <climits>

namespace detray {

/// Aborter that checks whether the track has exceeded its pathlimit
///
/// @tparam ID the actor id for this aborter ties its state instance
template<std::size_t ID>
struct pathlimit_aborter : actor<ID> {

    // Tag this actor
    using actor_type = pathlimit_aborter<ID>;

    // Pathlimit for a single propagation workflow
    struct state : actor<ID>::state {
        scalar _path_limit = std::numeric_limits<scalar>::max();

        /// Set the path limit to a scalar @param pl
        DETRAY_HOST_DEVICE
        inline void set_path_limit(const scalar pl) { _path_limit = pl; }

        /// @returns this states remaining path length.
        DETRAY_HOST_DEVICE
        inline scalar path_limit() const { return _path_limit; }
    }

    template <typename propagator_state_t>
    void operator()(actor_type::state & state, propagator_state_t& prop_state) {
        auto &step_state = prop_state._stepping;
        auto &nav_state = prop_state._navigation;

        // Nothing left to do. Propagation will exit successfully
        if (nav_state.is_complete()) {
            return;
        }

        state.path_limit() -= std::abs(_step_size);
        if (state.path_limit() <= 0.) {
            printf("Abort: Stepper above maximal path length!\n");
            // Stop navigation
            nav_state.abort();
        }

        step_state.template set_constraint<step::constraint::e_aborter>(state._path_limit);
    }
};

}  // namespace detray