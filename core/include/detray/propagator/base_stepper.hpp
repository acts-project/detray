/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/constrained_step.hpp"
#include "detray/tracks/tracks.hpp"

namespace detray {

namespace stepping {

enum class id {
    e_linear = 0,
    e_rk = 1,
};

}  // namespace stepping

/// Base stepper implementation
template <typename transform3_t, typename constraint_t, typename policy_t>
class base_stepper {

    public:
    using transform3_type = transform3_t;
    using free_track_parameters_type = free_track_parameters<transform3_t>;
    using bound_track_parameters_type = bound_track_parameters<transform3_t>;
    using matrix_operator = typename transform3_t::matrix_actor;
    using track_helper = detail::track_helper<matrix_operator>;

    using size_type = typename transform3_type::size_type;
    using scalar_type = typename transform3_type::scalar_type;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename transform3_type::matrix_actor::template matrix_type<ROWS,
                                                                     COLS>;
    /// Shorthand vector/matrix types related to bound track parameters.
    using bound_vector = matrix_type<e_bound_size, 1>;
    using bound_matrix = matrix_type<e_bound_size, e_bound_size>;
    /// Mapping from bound track parameters.
    using bound_to_free_matrix = matrix_type<e_free_size, e_bound_size>;
    // Shorthand vector/matrix types related to free track parameters.
    using free_vector = matrix_type<e_free_size, 1>;
    using free_matrix = matrix_type<e_free_size, e_free_size>;
    // Mapping from free track parameters.
    using free_to_bound_matrix = matrix_type<e_bound_size, e_free_size>;
    using free_to_path_matrix = matrix_type<1, e_free_size>;

    /** State struct holding the track
     *
     * It has to cast into a const track via the call
     * operation.
     */
    struct state {

        /// Sets track parameters.
        DETRAY_HOST_DEVICE
        state(const free_track_parameters_type &t) : _track(t) {}

        /// Sets track parameters from bound track parameter.
        template <typename detector_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type &bound_params,
            const detector_t &det)
            : _bound_params(bound_params) {

            const auto &surface_container = det.surfaces();
            const auto &trf_store = det.transform_store();
            const auto &mask_store = det.mask_store();
            const auto &surface =
                surface_container[bound_params.surface_link()];

            mask_store.template execute<
                typename parameter_resetter<transform3_t>::kernel>(
                surface.mask_type(), trf_store, surface, *this);
        }

        /// free track parameter
        free_track_parameters_type _track;

        /// jacobian transport matrix
        free_matrix _jac_transport =
            matrix_operator().template identity<e_free_size, e_free_size>();

        /// bound-to-free jacobian from departure surface
        bound_to_free_matrix _jac_to_global =
            matrix_operator().template zero<e_free_size, e_bound_size>();

        /// bound covariance
        bound_track_parameters_type _bound_params;

        /// @returns track parameters - const access
        DETRAY_HOST_DEVICE
        free_track_parameters_type &operator()() { return _track; }

        /// @returns track parameters.
        DETRAY_HOST_DEVICE
        const free_track_parameters_type &operator()() const { return _track; }

        step::direction _direction{step::direction::e_forward};

        // Stepping constraints
        constraint_t _constraint = {};

        // Navigation policy state
        typename policy_t::state _policy_state = {};

        /// Track path length
        scalar _path_length{0.};

        /// Track path length from the last surface. It will be reset to 0 when
        /// the track reaches a new surface
        scalar _s{0.};

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