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
template <typename track_t, template <typename...> class tuple_t = dtuple>
class line_stepper final : public base_stepper<track_t, tuple_t> {

    public:
    using base_type = base_stepper<track_t, tuple_t>;

    struct state : public base_type::state {
        state(track_t &t) : base_type::state(t) {}

        // Remaining path limit
        scalar _pl = std::numeric_limits<scalar>::max();

        /** @return the step and heartbeat given a step length s */
        tuple_t<scalar, bool> step(scalar s) {
            _pl = (s > _pl) ? s - _pl : _pl - s;
            const bool heartbeat = (_pl > 0.);
            return std::tie(std::min(s, _pl), heartbeat);
        }

        /** Set the path limit to a scalar @param l */
        void set_limit(scalar pl) { _pl = pl; }
    };

    /** Take a step, regulared by a constrained step
     *
     * @param s The state object that chaches
     * @param es The external step, e.g. from navigation
     *
     * @return returning the heartbeat, indicating if the stepping is alive
     */
    DETRAY_HOST_DEVICE
    bool step(state &s, scalar es) {
        const auto [sl, heartbeat] = s.step(es);
        s._track.set_pos(s._track.pos() + s._track.dir() * sl);
        return heartbeat;
    }
};

}  // namespace detray
