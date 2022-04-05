/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <climits>

#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Templated propagator class, using a stepper and a navigator object in
 *  succession.
 *
 * @tparam stepper_t for the transport
 * @tparam navigator_t for the navigation
 *
 **/
template <typename stepper_t, typename navigator_t, typename actor_chain_t>
struct propagator {

    stepper_t _stepper;
    navigator_t _navigator;

    // Register the actor types
    const actor_chain_t run_actors{};

    template <typename T>
    using vector_type = typename navigator_t::template vector_type<T>;

    /** Can not be default constructed */
    propagator() = delete;

    /** Only valid constructor with a
     * @param s stepper
     * @param n navigator
     * by move semantics
     **/
    DETRAY_HOST_DEVICE
    propagator(stepper_t &&s, navigator_t &&n)
        : _stepper(std::move(s)), _navigator(std::move(n)) {}

    /** Propagation that state aggregates a stepping and a navigation state */
    struct state {

        template <typename track_t>
        DETRAY_HOST_DEVICE state(
            track_t &t_in, typename actor_chain_t::state actor_states = {},
            vector_type<line_plane_intersection> &&candidates = {})
            : _stepping(t_in),
              _navigation(std::move(candidates)),
              _actor_states(std::move(actor_states)) {}

        typename stepper_t::state _stepping;
        typename navigator_t::state _navigation;
        typename actor_chain_t::state _actor_states;
    };

    /** Propagate method
     *
     * @tparam state_t is the propagation state type
     * @tparam inspector_t is the type of a propagation inspector
     *
     * @param p_state the state of a propagation
     * @param inspector the inspector
     *
     * @return propagation success.
     */
    template <typename state_t>
    DETRAY_HOST_DEVICE bool propagate(state_t &p_state) {

        auto &n_state = p_state._navigation;
        auto &s_state = p_state._stepping;
        auto &actor_states = p_state._actor_states;

        // initialize the navigation
        bool heartbeat = _navigator.init(n_state, s_state);

        // Run while there is a heartbeat
        while (heartbeat) {

            // Take the step
            heartbeat &= _stepper.step(s_state, n_state);

            // And check the status
            heartbeat &= _navigator.update(n_state, s_state);

            // Run all registered actors
            run_actors(actor_states, p_state);
        }

        return n_state.is_complete();
    }
};

}  // namespace detray