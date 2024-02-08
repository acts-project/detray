/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/macros.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagation_config.hpp"
#include "detray/tracks/tracks.hpp"

// System include(s).
#include <iomanip>

namespace detray {

/// Templated propagator class, using a stepper and a navigator object in
/// succession.
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
    using scalar_type = typename transform3_type::scalar_type;
    using free_track_parameters_type =
        typename stepper_t::free_track_parameters_type;
    using bound_track_parameters_type =
        typename stepper_t::bound_track_parameters_type;

    propagation::config<scalar_type> m_cfg;

    stepper_t m_stepper;
    navigator_t m_navigator;

    /// Register the actor types
    const actor_chain_t run_actors{};

    template <typename T>
    using vector_type = typename navigator_t::template vector_type<T>;

    /// Construct from a propagator configuration
    DETRAY_HOST_DEVICE
    propagator(propagation::config<scalar_type> cfg = {})
        : m_cfg{cfg}, m_stepper{}, m_navigator{} {}

    /// Propagation that state aggregates a stepping and a navigation state. It
    /// also keeps references to the actor states.
    struct state {

        using detector_type = typename navigator_t::detector_type;
        using navigator_state_type = typename navigator_t::state;
        using scalar_type = typename navigator_t::scalar_type;

        /// Construct the propagation state.
        ///
        /// @param t_in the track state to be propagated
        /// @param actor_states tuple that contains references to actor states
        /// @param candidates buffer for intersections in the navigator
        DETRAY_HOST_DEVICE state(
            const free_track_parameters_type &t_in, const detector_type &det,
            vector_type<intersection_type> &&candidates = {})
            : m_param_type(parameter_type::e_free),
              _stepping(t_in),
              _navigation(det, std::move(candidates)) {}

        template <typename field_t>
        DETRAY_HOST_DEVICE state(
            const free_track_parameters_type &t_in,
            const field_t &magnetic_field, const detector_type &det,
            vector_type<intersection_type> &&candidates = {})
            : m_param_type(parameter_type::e_free),
              _stepping(t_in, magnetic_field),
              _navigation(det, std::move(candidates)) {}

        /// Construct the propagation state with bound parameter
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type &param, const detector_type &det,
            vector_type<intersection_type> &&candidates = {})
            : m_param_type(parameter_type::e_bound),
              _stepping(param, det),
              _navigation(det, std::move(candidates)) {}

        /// Construct the propagation state with bound parameter
        template <typename field_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type &param,
            const field_t &magnetic_field, const detector_type &det,
            vector_type<intersection_type> &&candidates = {})
            : m_param_type(parameter_type::e_bound),
              _stepping(param, magnetic_field, det),
              _navigation(det, std::move(candidates)) {}

        DETRAY_HOST_DEVICE
        parameter_type param_type() const { return m_param_type; }

        DETRAY_HOST_DEVICE
        void set_param_type(const parameter_type t) { m_param_type = t; }

        // Is the propagation still alive?
        bool _heartbeat = false;
        // Starting type of track parametrization
        parameter_type m_param_type = parameter_type::e_free;

        typename stepper_t::state _stepping;
        typename navigator_t::state _navigation;

        bool do_debug = false;
#if defined(__NO_DEVICE__)
        std::stringstream debug_stream{};
#endif
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
        propagation._heartbeat =
            m_navigator.init(propagation, m_cfg.navigation);

        // Run all registered actors/aborters after init
        run_actors(actor_states, propagation);

        // Find next candidate
        propagation._heartbeat &=
            m_navigator.update(propagation, m_cfg.navigation);

        // Run while there is a heartbeat
        while (propagation._heartbeat) {

            // Take the step
            propagation._heartbeat &=
                m_stepper.step(propagation, m_cfg.stepping);

            // Find next candidate
            propagation._heartbeat &=
                m_navigator.update(propagation, m_cfg.navigation);

            // Run all registered actors/aborters after update
            run_actors(actor_states, propagation);

            // And check the status
            propagation._heartbeat &=
                m_navigator.update(propagation, m_cfg.navigation);

#if defined(__NO_DEVICE__)
            if (propagation.do_debug) {
                inspect(propagation);
            }
#endif
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
        propagation._heartbeat =
            m_navigator.init(propagation, m_cfg.navigation);

        // Run all registered actors/aborters after init
        run_actors(actor_states, propagation);

        // Find next candidate
        propagation._heartbeat &=
            m_navigator.update(propagation, m_cfg.navigation);

        while (propagation._heartbeat) {

            while (propagation._heartbeat) {

                // Take the step
                propagation._heartbeat &=
                    m_stepper.step(propagation, m_cfg.stepping);

                // Find next candidate
                propagation._heartbeat &=
                    m_navigator.update(propagation, m_cfg.navigation);

                // If the track is on a sensitive surface, break the loop to
                // synchornize the threads
                if (propagation._navigation.is_on_sensitive()) {
                    break;
                } else {
                    run_actors(actor_states, propagation);

                    // And check the status
                    propagation._heartbeat &=
                        m_navigator.update(propagation, m_cfg.navigation);
                }
            }

            // Synchornized actor
            if (propagation._heartbeat) {
                run_actors(actor_states, propagation);

                // And check the status
                propagation._heartbeat &=
                    m_navigator.update(propagation, m_cfg.navigation);

#if defined(__NO_DEVICE__)
                if (propagation.do_debug) {
                    inspect(propagation);
                }
#endif
            }
        }

        // Pass on the whether the propagation was successful
        return propagation._navigation.is_complete();
    }

    template <typename state_t>
    DETRAY_HOST void inspect(state_t &propagation) {
        const auto &navigation = propagation._navigation;
        const auto &stepping = propagation._stepping;

        propagation.debug_stream << std::left << std::setw(30);
        switch (navigation.status()) {
            case navigation::status::e_abort:
                propagation.debug_stream << "status: abort";
                break;
            case navigation::status::e_on_target:
                propagation.debug_stream << "status: e_on_target";
                break;
            case navigation::status::e_unknown:
                propagation.debug_stream << "status: unknowm";
                break;
            case navigation::status::e_towards_object:
                propagation.debug_stream << "status: towards_surface";
                break;
            case navigation::status::e_on_module:
                propagation.debug_stream << "status: on_module";
                break;
            case navigation::status::e_on_portal:
                propagation.debug_stream << "status: on_portal";
                break;
        };

        if (detail::is_invalid_value(navigation.volume())) {
            propagation.debug_stream << "volume: " << std::setw(10)
                                     << "invalid";
        } else {
            propagation.debug_stream << "volume: " << std::setw(10)
                                     << navigation.volume();
        }

        propagation.debug_stream << "surface: " << std::setw(14);
        if (navigation.is_on_portal() or navigation.is_on_module()) {
            propagation.debug_stream << navigation.barcode();
        } else {
            propagation.debug_stream << "undefined";
        }

        propagation.debug_stream << "step_size: " << std::setw(10)
                                 << stepping._step_size << std::endl;

        propagation.debug_stream << std::setw(10)
                                 << detail::ray<transform3_type>(stepping())
                                 << std::endl;
    }
};

}  // namespace detray
