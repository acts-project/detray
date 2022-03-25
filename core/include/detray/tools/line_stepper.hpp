/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray tools
#include "detray/tools/base_stepper.hpp"

// detray definitions
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Line stepper implementation */
template <typename track_t>
class line_stepper final : public base_stepper<track_t> {

    public:
    using base_type = base_stepper<track_t>;

    struct state : public base_type::state {
        state(track_t &t) : base_type::state(t) {}

        /// Accumulated path length
        scalar _path_length = 0.;

        /// Update the track state in a straight line.
        DETRAY_HOST_DEVICE
        inline void advance_track() {
            auto &track = this->_track;
            track.set_pos(track.pos() + track.dir() * this->_step_size);

            this->_path_length += this->_step_size;
        }
    };

    /** Take a step, regulared by a constrained step
     *
     * @param stepping The state object of a stepper
     * @param navigation The state object of a navigator
     * @param max_step_size Maximal distance for this step
     *
     * @return returning the heartbeat, indicating if the stepping is alive
     */
    template <typename navigation_state_t>
    DETRAY_HOST_DEVICE bool step(
        state &stepping, navigation_state_t &navigation,
        scalar max_step_size = std::numeric_limits<scalar>::max()) {

        // Distance to next surface as fixed step size
        scalar step_size = navigation();

        // Inform navigator
        // Not a severe change to track state expected
        if (step_size < max_step_size) {
            stepping.set_step_size(step_size);
            navigation.set_high_trust();
        }
        // Step size hit a constraint - the track state was probably changed a
        // lot
        else {
            stepping.set_step_size(max_step_size);
            // Re-evaluate all candidates
            navigation.set_fair_trust();
        }

        // Update and check path limit
        if (not stepping.check_path_limit()) {
            printf("Stepper: Above maximal path length!\n");
            // State is broken
            return navigation.abort();
        }

        // Update track state
        stepping.advance_track();

        return true;
    }
};

}  // namespace detray
