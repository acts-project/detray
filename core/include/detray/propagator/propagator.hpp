/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/tracks/tracks.hpp"

namespace detray {

/// Templated propagator class, using a stepper and a navigator object in
///  succession.
///
/// @tparam stepper_t for the transport
/// @tparam navigator_t for the navigation
template <typename stepper_t, typename navigator_t, typename actor_chain_t>
struct propagator {

    using stepper_type = stepper_t;
    using navigator_type = navigator_t;
    using intersection_type = typename navigator_type::intersection_type;
    using detector_type = typename navigator_type::detector_type;
    using actor_chain_type = actor_chain_t;
    using transform3_type = typename stepper_t::transform3_type;
    using free_track_parameters_type =
        typename stepper_t::free_track_parameters_type;
    using bound_track_parameters_type =
        typename stepper_t::bound_track_parameters_type;

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

        using detector_type = typename navigator_t::detector_type;
        using navigator_state_type = typename navigator_t::state;

        /// Construct the propagation state.
        ///
        /// @param t_in the track state to be propagated
        /// @param actor_states tuple that contains references to actor states
        /// @param candidates buffer for intersections in the navigator
        DETRAY_HOST_DEVICE state(
            const free_track_parameters_type &t_in, const detector_type &det,
            vector_type<intersection_type> &&candidates = {})
            : _stepping(t_in),
              _navigation(det, std::move(candidates)),
              m_param_type(parameter_type::e_free) {}

        template <typename field_t>
        DETRAY_HOST_DEVICE state(
            const free_track_parameters_type &t_in,
            const field_t &magnetic_field, const detector_type &det,
            vector_type<intersection_type> &&candidates = {})
            : _stepping(t_in, magnetic_field),
              _navigation(det, std::move(candidates)),
              m_param_type(parameter_type::e_free) {}

        /// Construct the propagation state with bound parameter
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type &param, const detector_type &det,
            vector_type<intersection_type> &&candidates = {})
            : _stepping(param, det),
              _navigation(det, std::move(candidates)),
              m_param_type(parameter_type::e_bound) {}

        /// Construct the propagation state with bound parameter
        template <typename field_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type &param,
            const field_t &magnetic_field, const detector_type &det,
            vector_type<intersection_type> &&candidates = {})
            : _stepping(param, magnetic_field, det),
              _navigation(det, std::move(candidates)),
              m_param_type(parameter_type::e_bound) {}

        DETRAY_HOST_DEVICE
        parameter_type param_type() const { return m_param_type; }

        DETRAY_HOST_DEVICE
        void set_param_type(const parameter_type t) { m_param_type = t; }

        // Is the propagation still alive?
        bool _heartbeat = false;

        typename stepper_t::state _stepping;
        typename navigator_t::state _navigation;
        parameter_type m_param_type = parameter_type::e_free;
    };

    /// Propagate method: Coordinates the calls of the stepper, navigator and
    /// all registered actors.
    ///
    /// @tparam state_t is the propagation state type
    /// @tparam actor_state_t is the actor state type
    ///
    /// @param propagation the state of a propagation flow
    /// @param actor_states the actor state
    ///
    /// @return propagation success.
    template <typename state_t, typename actor_states_t = actor_chain<>::state>
    DETRAY_HOST_DEVICE bool propagate(state_t &propagation,
                                      actor_states_t &&actor_states = {}) {

        // Initialize the navigation
        propagation._heartbeat = _navigator.init(propagation);

        // Run all registered actors/aborters after init
        run_actors(actor_states, propagation);

        // Run while there is a heartbeat
        while (propagation._heartbeat) {

            // Take the step
            propagation._heartbeat &= _stepper.step(propagation);

            // Find next candidate
            propagation._heartbeat &= _navigator.update(propagation);

            // Run all registered actors/aborters after update
            run_actors(actor_states, propagation);

            // And check the status
            propagation._heartbeat &= _navigator.update(propagation);
        }

        // Pass on the whether the propagation was successful
        return propagation._navigation.is_complete();
    }

    /// Propagate method with two while loops. In the CPU, propagate and
    /// propagate_sync() should be equivalent to each other. In the SIMT level
    /// (e.g. GPU), the instruction of threads in the same warp is synchornized
    /// after the internal while loop.
    ///
    /// @tparam state_t is the propagation state type
    /// @tparam actor_state_t is the actor state type
    ///
    /// @param propagation the state of a propagation flow
    /// @param actor_states the actor state
    ///
    /// @return propagation success.
    template <typename state_t, typename actor_states_t = actor_chain<>::state>
    DETRAY_HOST_DEVICE bool propagate_sync(state_t &propagation,
                                           actor_states_t &&actor_states = {}) {

        // Initialize the navigation
        propagation._heartbeat = _navigator.init(propagation);

        // Run all registered actors/aborters after init
        run_actors(actor_states, propagation);

        while (propagation._heartbeat) {

            while (propagation._heartbeat) {

                // Take the step
                propagation._heartbeat &= _stepper.step(propagation);

                // Find next candidate
                propagation._heartbeat &= _navigator.update(propagation);

                // If the track is on a sensitive surface, break the loop to
                // synchornize the threads
                if (propagation._navigation.is_on_sensitive()) {
                    break;
                } else {
                    run_actors(actor_states, propagation);

                    // And check the status
                    propagation._heartbeat &= _navigator.update(propagation);
                }
            }

            // Synchornized actor
            if (propagation._heartbeat) {
                run_actors(actor_states, propagation);

                // And check the status
                propagation._heartbeat &= _navigator.update(propagation);
            }
        }

        // Pass on the whether the propagation was successful
        return propagation._navigation.is_complete();
    }
};

}  // namespace detray
