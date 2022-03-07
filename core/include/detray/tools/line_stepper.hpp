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
#include <climits>

#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Line stepper implementation */
template <typename track_t>
class line_stepper final : public base_stepper<track_t> {

    public:
    using base_type = base_stepper<track_t>;

    struct state : public base_type::state {
        state(track_t &t) : base_type::state(t) {}

        // Remaining path limit
        scalar _path_limit = std::numeric_limits<scalar>::max();

        scalar _step_size = std::numeric_limits<scalar>::infinity();

        /** Set the path limit to a scalar @param l */
        void set_limit(scalar pl) { _path_limit = pl; }

        /// Set next step size
        DETRAY_HOST_DEVICE
        void set_step_size(const scalar /*step*/) {}
    };

    /** Take a step, regulared by a constrained step
     *
     * @param s The state object that chaches
     * @param es The external step, e.g. from navigation
     *
     * @return returning the heartbeat, indicating if the stepping is alive
     */
    template <typename navigation_state_t>
    DETRAY_HOST_DEVICE bool step(
        state &s, navigation_state_t &navigation,
        scalar max_step_size = std::numeric_limits<scalar>::max()) {
        scalar step_size = navigation();
        s._path_limit = (step_size > s._path_limit) ? step_size - s._path_limit
                                                    : s._path_limit - step_size;
        if (s._path_limit <= 0.) {
            return false;
        }
        s._step_size =
            std::min(step_size, std::min(s._path_limit, max_step_size));
        s._track.set_pos(s._track.pos() + s._track.dir() * s._step_size);
        navigation.set_high_trust();
        return true;
    }
};

}  // namespace detray
