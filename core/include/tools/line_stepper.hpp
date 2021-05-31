
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "utils/containers.hpp"

namespace detray
{

    /** Line stepper implementation */
    template <typename track_type>
    struct line_stepper {

        /** State struct holding the track 
         * 
         * It has to cast into a const track via the call
         * operat.
         */
        struct state {
           
            track_type& _track;

            scalar _s = 0.; //!< Next step 
            scalar _pl = std::numeric_limits<scalar>::max(); //!< Remaining path limit

            state() = delete;
            state(track_type& t) : _track(t) {}


            /** @return the step and heartbeat given a step length s */
            dtuple<scalar, bool> step(scalar s) {
                _pl = (s > _pl) ? s - _pl : _pl - s;
                const bool heartbeat = (_pl > 0.);
                return std::tie(std::min(s,_pl), heartbeat);
            }

            /** Set the path limit to a scalar @param l */
            void set_limit(scalar pl) {
                _pl = pl;
            }

            /** Call operator casts it itno a const track referenc on @return */
            const track_type& operator()() const { return _track; }

        };

        /** Take a step, regulared by a constrained step
         * 
         * @param s The state object that chaches
         * @param es The external step, e.g. from navigation 
         * 
         * @return returning the heartbeat, indicating if the stepping is alive
         */
        bool step(state& s, scalar es) const {
            const auto [ sl, heartbeat ] = s.step(es);
            s._track.pos = s._track.pos + s._track.dir * sl;
            return heartbeat;
        }

    };

}
