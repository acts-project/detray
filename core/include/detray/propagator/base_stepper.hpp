/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/constrained_step.hpp"
#include "detray/propagator/detail/covariance_engine.hpp"
#include "detray/propagator/track.hpp"

namespace detray {

/// Base stepper implementation
///
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename track_t, typename constraint_t, typename policy_t>
class base_stepper {

    public:
    // covariance types for bound state
    using covariance_type = typename bound_track_parameters::covariance_type;
    using jacobian_type = typename bound_track_parameters::jacobian_type;
    using matrix_operator = standard_matrix_operator<scalar>;
    using covariance_engine = detail::covariance_engine<scalar>;
    using vector_engine = covariance_engine::vector_engine;
    using transform3 = typename covariance_engine::transform3;

    /** State struct holding the track
     *
     * It has to cast into a const track via the call
     * operation.
     */
    struct state {

        /// Sets track parameters.
        DETRAY_HOST_DEVICE
        state(const track_t &t) : _track(t) {}

        /// Sets track parameters from bound track parameter.
        DETRAY_HOST_DEVICE
        state(const bound_track_parameters &bound_params,
              const transform3 &trf3) {
            // Set the free vector
            _track.set_vector(vector_engine().bound_to_free_vector(
                trf3, bound_params.vector()));

            // Set the bound covariance
            _bound_covariance = bound_params.covariance();

            // Reset the jacobians
            covariance_engine().reinitialize_jacobians(
                trf3, bound_params.vector(), _jac_to_global, _jac_transport,
                _derivative);
        }

        /// free track parameter
        track_t _track;

        /// jacobian transport matrix
        free_matrix _jac_transport =
            matrix_operator().template identity<e_free_size, e_free_size>();

        /// The free parameter derivative defined at destination surface
        free_vector _derivative =
            matrix_operator().template zero<e_free_size, 1>();

        /// bound-to-free jacobian from departure surface
        bound_to_free_matrix _jac_to_global =
            matrix_operator().template zero<e_free_size, e_bound_size>();

        /// bound covariance
        bound_matrix _bound_covariance =
            matrix_operator().template zero<e_bound_size, e_bound_size>();

        /// @returns track parameters - const access
        DETRAY_HOST_DEVICE
        track_t &operator()() { return _track; }

        /// @returns track parameters.
        DETRAY_HOST_DEVICE
        const track_t &operator()() const { return _track; }

        step::direction _direction{step::direction::e_forward};

        // Stepping constraints
        constraint_t _constraint = {};

        // Navigation policy state
        typename policy_t::state _policy_state = {};

        /// Track path length
        scalar _path_length{0.};

        /// Current step size
        scalar _step_size{0.};

        /// TODO: Use options?
        /// hypothetical mass of particle (assume pion by default)
        /// scalar _mass = 139.57018 * unit_constants::MeV;

        /// Set new step constraint
        template <step::constraint type = step::constraint::e_actor>
        DETRAY_HOST_DEVICE inline void set_constraint(scalar step_size) {
            _constraint.template set<type>(step_size);
        }

        /// Set new navigation direction
        DETRAY_HOST_DEVICE inline void set_direction(step::direction dir) {
            _direction = dir;
        }

        /// @returns access to this states step constraints
        DETRAY_HOST_DEVICE
        inline const constraint_t &constraints() const { return _constraint; }

        /// @returns access to this states step constraints
        DETRAY_HOST_DEVICE
        inline typename policy_t::state &policy_state() {
            return _policy_state;
        }

        /// @returns the navigation direction
        DETRAY_HOST_DEVICE
        inline step::direction direction() const { return _direction; }

        /// Remove [all] constraints
        template <step::constraint type = step::constraint::e_actor>
        DETRAY_HOST_DEVICE inline void release_step() {
            _constraint.template release<type>();
        }

        /// Set next step size
        DETRAY_HOST_DEVICE
        inline void set_step_size(const scalar step) { _step_size = step; }

        /// @returns the current step size of this state.
        DETRAY_HOST_DEVICE
        inline scalar step_size() const { return _step_size; }

        /// @returns this states remaining path length.
        DETRAY_HOST_DEVICE
        inline scalar path_length() const { return _path_length; }
    };
};

}  // namespace detray