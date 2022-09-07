/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/base_stepper.hpp"

namespace detray {

struct surface_targeter : actor {

    struct state {

        scalar _path;
        dindex _target_surface_index = dindex_invalid;
    };

    /// Enforces thepath limit on a stepper state
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

        scalar residual = actor_state._path - stepping.path_length();
        stepping.set_constraint(residual);

        typename propagator_state_t::navigator_state_type::intersection_t is;
        is.index = actor_state._target_surface_index;
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