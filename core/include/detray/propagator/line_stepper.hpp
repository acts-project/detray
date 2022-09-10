/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray definitions
#include "detray/definitions/qualifiers.hpp"

// detray tools
#include "detray/propagator/base_stepper.hpp"
#include "detray/propagator/navigation_policies.hpp"

// system includes
#include <cmath>

namespace detray {

/// Straight line stepper implementation
template <typename transform3_t, typename constraint_t = unconstrained_step,
          typename policy_t = stepper_default_policy>
class line_stepper final
    : public base_stepper<transform3_t, constraint_t, policy_t> {

    public:
    using base_type = base_stepper<transform3_t, constraint_t, policy_t>;
    using policy_type = policy_t;
    using free_track_parameters_type =
        typename base_type::free_track_parameters_type;
    using bound_track_parameters_type =
        typename base_type::bound_track_parameters_type;
    using transform3_type = transform3_t;

    struct state : public base_type::state {

        static constexpr const stepping::id id = stepping::id::e_linear;

        DETRAY_HOST_DEVICE
        state(const free_track_parameters_type &t) : base_type::state(t) {}

        template <typename detector_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type &bound_params,
            const detector_t &det)
            : base_type::template state(bound_params, det) {}

        /// Update the track state in a straight line.
        DETRAY_HOST_DEVICE
        inline void advance_track() {
            auto &track = this->_track;
            track.set_pos(track.pos() + track.dir() * this->_step_size);

            this->_path_length += this->_step_size;
        }
    };

    /// Take a step, regulared by a constrained step
    ///
    /// @param stepping The state object of a stepper
    /// @param navigation The state object of a navigator
    /// @param max_step_size Maximal distance for this step
    ///
    /// @return returning the heartbeat, indicating if the stepping is alive
    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE bool step(propagation_state_t &propagation) {
        // Get stepper state
        state &stepping = propagation._stepping;
        // Distance to next surface as fixed step size
        scalar step_size = propagation._navigation();

        // Update navigation direction
        const step::direction dir = step_size > 0 ? step::direction::e_forward
                                                  : step::direction::e_backward;
        stepping.set_direction(dir);

        // Check constraints
        if (std::abs(step_size) >
            std::abs(
                stepping.constraints().template size<>(stepping.direction()))) {
            stepping.set_step_size(
                stepping.constraints().template size<>(stepping.direction()));
        } else {
            stepping.set_step_size(step_size);
        }

        // Update track state
        stepping.advance_track();

        // Call navigation update policy
        policy_t{}(stepping.policy_state(), propagation);

        return true;
    }
};

}  // namespace detray
