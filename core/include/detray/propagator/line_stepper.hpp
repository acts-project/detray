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
///
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename track_t, typename constraint_t = unconstrained_step,
          typename policy_t = stepper_default_policy>
class line_stepper final
    : public base_stepper<track_t, constraint_t, policy_t> {

    public:
    using base_type = base_stepper<track_t, constraint_t, policy_t>;
    using policy_type = policy_t;

    struct state : public base_type::state {

        using field_type = int;

        state(){};

        DETRAY_HOST_DEVICE
        state(const track_t &t) : base_type::state(t) {}

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
