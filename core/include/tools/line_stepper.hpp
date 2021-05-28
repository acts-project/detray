
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray
{

    /** Line stepper implementation */
    template <typename track_type>
    struct line_stepper {

        /** State struct holding the track */
        struct state {
            state() = delete;
            state(track_type& t) : _track(t) {}

            track_type& _track;

            scalar _s = 0.;
        };

        /** Take a step, regulared by a constrained step
         * 
         * @param s The state object that chaches
         * @param cs The constraint step 
         * 
         * @return a boolean for abort 
         */
        bool step(state& s){
            s._track.pos = s._track.pos + s._track.dir * s._s;
            return false;
        }

    };

}
