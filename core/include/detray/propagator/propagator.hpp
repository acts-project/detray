/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <climits>
#include <iostream>

#include "detray/definitions/qualifiers.hpp"

namespace detray {

/// Templated propagator class, using a stepper and a navigator object in
///  succession.
///
/// @tparam stepper_t for the transport
/// @tparam navigator_t for the navigation
template <typename stepper_t, typename navigator_t, typename actor_chain_t>
struct propagator {

    stepper_t _stepper;
    navigator_t _navigator;

    /// Register the actor types
    const actor_chain_t run_actors{};

    template <typename T>
    using vector_type = typename navigator_t::template vector_type<T>;

    /// Cannot be default constructed
    propagator() = delete;

    /// Only valid constructor with a
    /// @param s stepper
    /// @param n navigator
    /// by move semantics
    DETRAY_HOST_DEVICE
    propagator(stepper_t &&s, navigator_t &&n)
        : _stepper(std::move(s)), _navigator(std::move(n)) {}

    /// Propagation that state aggregates a stepping and a navigation state. It
    /// also keeps references to the actor states.
    struct state {

        /// Construct the propagation state.
        ///
        /// @param t_in the track state to be propagated
        /// @param actor_states tuple that contains references to actor states
        /// @param candidates buffer for intersections in the navigator
        template <typename track_t>
        DETRAY_HOST_DEVICE state(
            track_t &t_in, typename actor_chain_t::state actor_states = {},
            vector_type<line_plane_intersection> &&candidates = {})
            : _stepping(t_in),
              _navigation(std::move(candidates)),
              _actor_states(actor_states) {}

        // Is the propagation still alive?
        bool _heartbeat = false;

        typename stepper_t::state _stepping;
        typename navigator_t::state _navigation;
        typename actor_chain_t::state _actor_states;
    };

    /// Propagate method: Coordinates the calls of the stepper, navigator and
    /// all registered actors.
    ///
    /// @tparam state_t is the propagation state type
    ///
    /// @param propagation the state of a propagation flow
    ///
    /// @return propagation success.
    template <typename state_t>
    DETRAY_HOST_DEVICE bool propagate(state_t &propagation) {

        // initialize the navigation
        propagation._heartbeat = _navigator.init(propagation);

        // Run while there is a heartbeat
        while (propagation._heartbeat) {
            // std::cout << propagation._navigation() << std::endl;

            /*if (propagation._navigation() == -0) {
                std::cout << propagation._navigation.inspector().to_string() <<
            std::endl; return false;
            }*/
            // Take the step
            propagation._heartbeat &= _stepper.step(propagation);

            // And check the status
            propagation._heartbeat &= _navigator.update(propagation);

            // Run all registered actors
            run_actors(propagation._actor_states, propagation);
        }

        // Pass on the whether the propagation was successful
        return propagation._navigation.is_complete();
    }
};

}  // namespace detray