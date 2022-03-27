/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray definitions
#include "detray/definitions/qualifiers.hpp"

// detray tools
#include "detray/tools/constrained_step.hpp"
#include "detray/tools/track.hpp"

namespace detray {

namespace step {

/// Direction in which the integration is performed
enum direction : int {
    e_forward = 1,
    e_unknown = std::numeric_limits<int>::max(),
    e_backward = -1,
};

}  // namespace step

/// Base stepper implementation
///
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename track_t, typename constraint_t>
class base_stepper {

    public:
    // covariance types for bound state
    using covariance_type = typename bound_track_parameters::covariance_type;
    using jacobian_type = typename bound_track_parameters::jacobian_type;

    /** State struct holding the track
     *
     * It has to cast into a const track via the call
     * operation.
     */
    struct state {

        state() = delete;

        /// Sets track parameters.
        DETRAY_HOST_DEVICE
        state(track_t &t) : _track(t) {}

        /// free track parameter
        track_t &_track;

        /// jacobian
        bound_matrix _jacobian;

        /// jacobian transport matrix
        free_matrix _jac_transport;

        /// jacobian transformation
        bound_to_free_matrix _jac_to_global;

        /// covariance matrix on surface
        bound_matrix _cov;

        /// The propagation derivative
        free_vector _derivative;

        /// @returns track parameters - const access
        DETRAY_HOST_DEVICE
        track_t &operator()() { return _track; }

        /// @returns track parameters.
        DETRAY_HOST_DEVICE
        const track_t &operator()() const { return _track; }

        step::direction _direction = step::direction::e_forward;

        // Stepping constraints
        constraint_t _constraint = {};

        /// Remaining path length
        scalar _path_limit = std::numeric_limits<scalar>::max();

        /// Current step size
        scalar _step_size = std::numeric_limits<scalar>::infinity();

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

        /// @returns the navigation direction
        DETRAY_HOST_DEVICE
        inline step::direction direction() const { return _direction; }

        /// Remove [all] constraints
        template <step::constraint type = step::constraint::e_actor>
        DETRAY_HOST_DEVICE inline void release_step() {
            _constraint.template release<type>();
        }

        /// Set the path limit to a scalar @param pl
        DETRAY_HOST_DEVICE
        inline void set_path_limit(const scalar pl) { _path_limit = pl; }

        /// @returns this states remaining path length.
        DETRAY_HOST_DEVICE
        inline scalar dist_to_path_limit() const { return _path_limit; }

        /// Update and check the path limit against a new @param step size.
        DETRAY_HOST_DEVICE
        inline bool check_path_limit() const {
            _path_limit -= _step_size;
            if (_path_limit <= 0.) {
                return false;
            }
            return true;
        }

        /// Set next step size
        DETRAY_HOST_DEVICE
        inline void set_step_size(const scalar step) { _step_size = step; }

        /// @returns the current step size of this state.
        DETRAY_HOST_DEVICE
        inline scalar step_size() const { return _step_size; }
    };
};

}  // namespace detray