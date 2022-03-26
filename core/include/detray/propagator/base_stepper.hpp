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
#include "detray/propagator/track.hpp"

namespace detray {

/** abstract stepper implementation */
template <typename track_t>
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

        enum navigation_direction : int {
            e_forward = 1,
            e_backward = -1,
        };

        navigation_direction _nav_dir = e_forward;

        /// Remaining path length
        scalar _path_limit = std::numeric_limits<scalar>::max();

        /// Current step size
        scalar _step_size = std::numeric_limits<scalar>::infinity();

        /// @returns this states remaining path length.
        DETRAY_HOST_DEVICE
        inline scalar dist_to_path_limit() const { return _path_limit; }

        /// Set the path limit to a scalar @param pl
        DETRAY_HOST_DEVICE
        inline void set_path_limit(const scalar pl) { _path_limit = pl; }

        /// Update and check the path limit against a new @param step size.
        DETRAY_HOST_DEVICE
        inline bool check_path_limit() {
            _path_limit -= _step_size;
            if (_path_limit <= 0.) {
                return false;
            }
            return true;
        }

        /// @returns the current step size of this state.
        DETRAY_HOST_DEVICE
        inline scalar step_size() const { return _step_size; }

        /// Set next step size
        DETRAY_HOST_DEVICE
        inline void set_step_size(const scalar step) { _step_size = step; }

        /// Set next step size and release constrained stepping
        DETRAY_HOST_DEVICE
        void release_step_size() {}
    };
};

}  // namespace detray