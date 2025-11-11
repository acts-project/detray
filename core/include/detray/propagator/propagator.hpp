/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/macros.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/propagator/concepts.hpp"
#include "detray/propagator/detail/noise_estimation.hpp"
#include "detray/propagator/propagation_config.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/log.hpp"

// System include(s).
#include <iomanip>

namespace detray {

/// Templated propagator class, using a stepper and a navigator object in
/// succession.
///
/// @tparam stepper_t for the transport
/// @tparam navigator_t for the navigation
template <typename stepper_t, typename navigator_t,
          concepts::actor_chain actor_chain_t>
struct propagator {

    using stepper_type = stepper_t;
    using navigator_type = navigator_t;
    using actor_chain_type = actor_chain_t;

    using detector_type = typename navigator_type::detector_type;
    using algebra_type = typename detector_type::algebra_type;
    using scalar_type = dscalar<algebra_type>;
    using intersection_type = typename navigator_type::intersection_type;
    using free_track_parameters_type =
        typename stepper_t::free_track_parameters_type;
    using bound_track_parameters_type =
        typename stepper_t::bound_track_parameters_type;

    propagation::config m_cfg;

    stepper_t m_stepper;
    navigator_t m_navigator;

    /// Register the actor types
    const actor_chain_type run_actors{};

    /// Construct from a propagator configuration
    DETRAY_HOST_DEVICE
    explicit constexpr propagator(const propagation::config &cfg)
        : m_cfg{cfg} {}

    /// Propagation that state aggregates a stepping and a navigation state. It
    /// also keeps references to the actor states.
    struct state {

        using detector_type = typename navigator_t::detector_type;
        using algebra_type = typename detector_type::algebra_type;
        using context_type = typename detector_type::geometry_context;
        using navigator_state_type = typename navigator_t::state;
        using actor_chain_type = actor_chain_t;
        using scalar_type = typename detector_type::scalar_type;

        /// Construct the propagation state with free parameter
        DETRAY_HOST_DEVICE state(const free_track_parameters_type &free_params,
                                 const detector_type &det,
                                 const context_type &ctx)
            : _stepping(free_params), _navigation(det), _context(ctx) {}

        /// Construct the propagation state with free parameter
        template <typename field_t>
        DETRAY_HOST_DEVICE state(const free_track_parameters_type &free_params,
                                 const field_t &magnetic_field,
                                 const detector_type &det,
                                 const context_type &ctx = {})
            : _stepping(free_params, magnetic_field),
              _navigation(det),
              _context(ctx) {}

        /// Construct the propagation state from the navigator state view
        DETRAY_HOST_DEVICE state(
            const free_track_parameters_type &free_params,
            const detector_type &det,
            typename navigator_type::state::view_type nav_view,
            const context_type &ctx = {})
            : _stepping(free_params),
              _navigation(det, nav_view),
              _context(ctx) {}

        /// Construct the propagation state from the navigator state view
        template <typename field_t>
        DETRAY_HOST_DEVICE state(
            const free_track_parameters_type &free_params,
            const field_t &magnetic_field, const detector_type &det,
            typename navigator_type::state::view_type nav_view,
            const context_type &ctx = {})
            : _stepping(free_params, magnetic_field),
              _navigation(det, nav_view),
              _context(ctx) {}

        /// Construct the propagation state with bound parameter
        DETRAY_HOST_DEVICE state(const bound_track_parameters_type &param,
                                 const detector_type &det,
                                 const context_type &ctx = {})
            : _stepping(param, det, ctx), _navigation(det), _context(ctx) {
            _navigation.set_volume(param.surface_link().volume());
        }

        /// Construct the propagation state with bound parameter
        template <typename field_t>
        DETRAY_HOST_DEVICE state(const bound_track_parameters_type &param,
                                 const field_t &magnetic_field,
                                 const detector_type &det,
                                 const context_type &ctx = {})
            : _stepping(param, magnetic_field, det, ctx),
              _navigation(det),
              _context(ctx) {
            _navigation.set_volume(param.surface_link().volume());
        }

        /// Construct the propagation state with bound parameter and navigator
        /// state view
        template <typename field_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type &param,
            const field_t &magnetic_field, const detector_type &det,
            typename navigator_type::state::view_type nav_view,
            const context_type &ctx = {})
            : _stepping(param, magnetic_field, det, ctx),
              _navigation(det, nav_view),
              _context(ctx) {
            _navigation.set_volume(param.surface_link().volume());
        }

        /// Set the particle hypothesis
        DETRAY_HOST_DEVICE
        void set_particle(const pdg_particle<scalar_type> &ptc) {
            _stepping.set_particle(ptc);
        }

        /// @returns the propagation heartbeat
        DETRAY_HOST_DEVICE
        bool is_alive() const { return _heartbeat; }

        // Is the propagation still alive?
        bool _heartbeat = false;

        typename stepper_t::state _stepping;
        typename navigator_t::state _navigation;
        context_type _context;

        bool do_debug = false;
    };

    /// Propagate method finale: Return whether or not the propagation
    /// completed succesfully.
    ///
    /// @param propagation the state of a propagation flow
    ///
    /// @return propagation success.
    DETRAY_HOST_DEVICE bool finished(const state &propagation) const {
        return propagation._navigation.finished();
    }

    /// @returns true if the @param propagation is suspended
    DETRAY_HOST_DEVICE
    inline auto is_paused(const state &propagation) const -> bool {
        return !propagation.is_alive() && propagation._navigation.is_alive();
    }

    /// Revive the propagation
    DETRAY_HOST_DEVICE
    inline void resume(state &propagation) const {
        assert(propagation._navigation.is_alive());
        propagation._heartbeat = true;
    }

    /// Propagate method: Coordinates the calls of the stepper, navigator and
    /// all registered actors.
    ///
    /// @param propagation the state of a propagation flow
    /// @param actor_state_refs tuple containing refences to the actor states
    ///
    /// @return propagation success.
    template <typename actor_states_t>
        requires concepts::is_state_of<actor_states_t, actor_chain_type>
    DETRAY_HOST_DEVICE bool propagate(
        state &propagation,
        actor_states_t actor_state_refs = dtuple<>{}) const {

        auto &navigation = propagation._navigation;
        auto &stepping = propagation._stepping;
        auto &context = propagation._context;
        const auto &track = stepping();
        assert(!track.is_invalid());

        DETRAY_VERBOSE_HOST("Starting propagation for track:\n" << track);

        // Open the navigation area according to uncertainties in initital track
        // params
        if (m_cfg.navigation.estimate_scattering_noise &&
            !stepping.bound_params().is_invalid()) {
            detail::estimate_external_mask_tolerance(
                stepping.bound_params(), propagation,
                static_cast<scalar_type>(m_cfg.navigation.n_scattering_stddev));
        }

        // Initialize the navigation
        m_navigator.init(track, navigation, m_cfg.navigation, context);
        propagation._heartbeat = navigation.is_alive();

        bool is_init = true;

        // Run while there is a heartbeat. In order to help the compiler
        // optimize this, and in order to make the code more GPU-friendly,
        // this code is run as a flat loop, but this loop has a defined
        // structure. Indeed, the structure is always to run either the actors
        // or the stepper (in alternating order) followed by the navigation
        // update.
        //
        // A = actors
        // N = navigation update
        // S = propagation step
        //
        // ANSNANSNANSNANSNANSNANS...
        scalar_type path_length{0.f};
        unsigned int stall_counter{0u};
        for (unsigned int i = 0; i % 2 == 0 || propagation.is_alive(); ++i) {

            if (i % 2 == 0) {
                DETRAY_VERBOSE_HOST_DEVICE("Propagation step: %d", i / 2);
                DETRAY_VERBOSE_HOST_DEVICE("Path length: %f mm",
                                           stepping.path_length());

                // Run all registered actors/aborters
                run_actors(actor_state_refs, propagation);

                // Don't run another navigation update, if already exited
                if (!propagation.is_alive()) {
                    continue;
                }

                path_length = stepping.path_length();

                assert(!track.is_invalid());
            } else {
                assert(!track.is_invalid());

                // Set access to the volume material for the stepper
                auto vol = navigation.current_volume();
                const material<scalar_type> *vol_mat_ptr =
                    vol.has_material() ? vol.material_parameters(track.pos())
                                       : nullptr;

                // Break automatic step size scaling by the stepper when a
                // surface was reached and whenever the navigation is
                // (re-)initialized
                const bool reset_stepsize{navigation.is_on_surface() ||
                                          is_init};

                // Take the step
                propagation._heartbeat &=
                    m_stepper.step(navigation(), stepping, m_cfg.stepping,
                                   reset_stepsize, vol_mat_ptr);

                // Reduce navigation trust level according to stepper update
                typename stepper_t::policy_type{}(stepping.policy_state(),
                                                  propagation);

                if (i > 0) {
                    is_init = false;
                }

                // Check if the propagation makes progress
                if (math::fabs(stepping.path_length()) <=
                    math::fabs(path_length +
                               m_cfg.navigation.intersection.path_tolerance)) {
                    if (stall_counter >= 10u) {
                        propagation._heartbeat = false;
                        navigation.abort("Propagation stalled");
                        DETRAY_ERROR_HOST("Propagation stalled");
                    } else if (stall_counter > 2u) {
                        // Print a warning if the propagation starts stalling
                        // (no overlap)
                        DETRAY_WARN_HOST_DEVICE(
                            "Propagation is stalling (counter %d)",
                            stall_counter);
                        DETRAY_WARN_HOST(print(propagation));
                        DETRAY_WARN_HOST("-> Track: " << stepping());
                    }
                    DETRAY_DEBUG_HOST_DEVICE("Step stalled. Counter %d",
                                             stall_counter);
                    stall_counter++;
                } else {
                    stall_counter = 0u;
                }
            }

            // Find next candidate
            is_init |= m_navigator.update(track, navigation, m_cfg.navigation,
                                          context, i % 2 == 1 || i == 0);

            propagation._heartbeat &= navigation.is_alive();

            if (i % 2 == 0 && i > 0 && propagation.do_debug) {
                DETRAY_VERBOSE_HOST(print(propagation));
            }
        }

        // Pass on the whether the propagation was successful
        DETRAY_VERBOSE_HOST("Finished propagation for track:\n" << track);
        if (finished(propagation)) {
            DETRAY_VERBOSE_HOST_DEVICE("Status: SUCCESS");
        } else if (is_paused(propagation)) {
            DETRAY_VERBOSE_HOST_DEVICE("Status: PAUSED");
        } else {
            DETRAY_VERBOSE_HOST_DEVICE("Status: ABORT");
        }

        return finished(propagation) || is_paused(propagation);
    }

    /// Overload for emtpy actor chain
    DETRAY_HOST_DEVICE bool propagate(state &propagation) {
        // Will not be used
        actor_chain<>::state empty_state{};
        // Run propagation
        return propagate(propagation, empty_state);
    }

    template <typename state_t>
    DETRAY_HOST std::string print(state_t &propagation) const {
        const auto &navigation = propagation._navigation;
        const auto &stepping = propagation._stepping;

        std::stringstream debug_stream{};
        debug_stream << std::left << std::setw(10);
        debug_stream << "status: " << navigation.status() << std::endl;

        debug_stream << "volume: " << std::setw(10);
        if (detail::is_invalid_value(navigation.volume())) {
            debug_stream << "invalid";
        } else {
            debug_stream << navigation.volume();
        }
        debug_stream << std::endl;

        debug_stream << "navigation:" << std::endl;
        if (navigation.is_on_surface()) {
            debug_stream << std::setw(10)
                         << " ->on surface: " << navigation.barcode();
        } else {
            debug_stream << std::setw(10) << " ->target: "
                         << navigation.target().sf_desc.barcode();
        }
        debug_stream << std::endl;
        debug_stream << " ->path: " << navigation() << "mm" << std::endl;

        debug_stream << "stepping:" << std::endl;
        debug_stream << " -> step size: " << std::setw(10)
                     << stepping.step_size() << "mm" << std::endl;
        debug_stream << " ->" << detail::ray<algebra_type>(stepping());

        return debug_stream.str();
    }
};

}  // namespace detray
