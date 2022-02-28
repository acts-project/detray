/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray definitions
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"

// detray tools
#include "detray/tools/track.hpp"

namespace detray {

/** abstract stepper implementation */
template <typename track_t, template <typename...> class tuple_t = dtuple>
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

        DETRAY_HOST_DEVICE
        state(track_t &t) : _track(t) {}

        // free track parameter
        track_t &_track;

        // jacobian
        bound_matrix _jacobian;

        // jacobian transport matrix
        free_matrix _jac_transport;

        // jacobian transformation
        bound_to_free_matrix _jac_to_global;

        // covariance matrix on surface
        bound_matrix _cov;

        // The propagation derivative
        free_vector _derivative;

        DETRAY_HOST_DEVICE
        track_t &operator()() { return _track; }

        DETRAY_HOST_DEVICE
        const track_t &operator()() const { return _track; }

        enum navigation_direction : int {
            e_forward = 1,
            e_backward = -1,
        };

        navigation_direction _nav_dir = e_forward;

        DETRAY_HOST_DEVICE
        void release_step_size(const scalar /*step*/) {}
    };
};

}  // namespace detray