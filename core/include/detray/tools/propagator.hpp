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

    /** Can not be default constructed */
    propagator() = delete;
    /** Only valid constructor with a
     * @param s stepper
     * @param n navigator
     * by move semantics
     **/
    propagator(stepper_t &&s, navigator_t &&n)
        : _stepper(std::move(s)), _navigator(std::move(n)) {}

    /** Propagate method
     *
     * @tparam track_t is the type of the track
     *
     * @param t_in the track at input
     *
     * @return a track at output
     */
    template <typename track_t, typename propagator_inspector_t>
    track_t propagate(const track_t &t_in,
                      propagator_inspector_t &propagator_inspector) {

        track_t t_out(t_in);
        typename stepper_t::state s_state(t_out);
        typename navigator_t::state n_state;
        // For now, always start at zero
        n_state.set_volume(0u);

        bool heartbeat = _navigator.status(n_state, s_state);

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
        return t_out;
    }
};

}  // namespace detray
