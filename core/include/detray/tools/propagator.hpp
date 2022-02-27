/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** A sterily track inspector instance */
struct void_propagator_inspector {

    /** void operator **/
    template <typename... args>
    DETRAY_HOST_DEVICE void operator()(const args &... /*ignored*/) {
        return;
    }
};

/** Tempalted propagator class, using a
 *
 * @tparam stepper_t for the transport
 * @tparam navigator_t for the navigation
 *
 **/
template <typename stepper_t, typename navigator_t>
struct propagator {

    stepper_t _stepper;
    navigator_t _navigator;

    template <typename T>
    using vector_type = typename navigator_t::template vector_type<T>;

    /** Can not be default constructed */
    propagator() = delete;
    /** Only valid constructor with a
     * @param s stepper
     * @param n navigator
     * by move semantics
     **/
    propagator(stepper_t &&s, navigator_t &&n)
        : _stepper(std::move(s)), _navigator(std::move(n)) {}

    struct state {

        template <typename track_t>
        state(track_t &t_in, vector_type<intersection> candidates = {})
            : _stepping(t_in), _navigation(candidates) {}

        typename stepper_t::state _stepping;
        typename navigator_t::state _navigation;
    };

    /** Propagate method
     *
     * @tparam track_t is the type of the track
     *
     * @param t_in the track at input
     *
     * @return a track at output
     */
    template <typename state_t, typename propagator_inspector_t>
    DETRAY_HOST_DEVICE void propagate(
        state_t &p_state, propagator_inspector_t &propagator_inspector) {

        auto &n_state = p_state._navigation;
        auto &s_state = p_state._stepping;

        // For now, always start at zero
        n_state.set_volume(0u);

        // bool heartbeat = _navigator.status(n_state, s_state);
        bool heartbeat = true;

        // Run while there is a heartbeat
        while (heartbeat) {

            // (Re-)target
            heartbeat &= _navigator.target(n_state, s_state);

            // Take the step
            heartbeat &= _stepper.step(s_state, n_state());

            // And check the status
            heartbeat &= _navigator.status(n_state, s_state);

            propagator_inspector(n_state, s_state);
        }
    }
};

}  // namespace detray
