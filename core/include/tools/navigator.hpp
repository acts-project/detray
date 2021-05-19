
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/intersection.hpp"
#include "core/track.hpp"
#include "utils/indexing.hpp"

namespace detray
{

    /** A void inpector that does nothing.
     * 
     * Inspectors can be plugged in to understand the
     * the current navigation state information.
     * 
     */
    struct void_inspector {
        template <typename state_type>
        void operator()(const state_type& /*ignored*/){}
    };


    /** The navigator struct. 
     * 
     * It follows the Acts::Navigator sturcture of sequence of
     * - status()
     * - target()
     * calls.
     * 
     * @tparam detector_tupe is the type of the detector
     * @tparam inspector_type is a validation inspector 
     */
    template <typename detector_type, typename inspector_type = void_inspector>
    struct navigator
    {

        using surface = typename detector_type::surface;
        using surface_link = typename detector_type::surface_link;

        using portal = typename detector_type::portal;
        using portal_links = typename detector_type::portal_links;

        using context = typename detector_type::context;

        /** Simple navigation status struct */
        enum navigation_status : int
        {
            e_on_target = -3,
            e_abort = -2,
            e_unknown = -1,
            e_towards_surface = 0,
            e_on_surface = 1,
            e_towards_portal = 2,
            e_on_portal = 3,
        };

        /** A nested navigation kernel struct which can be used for surfaces, portals,
         *  volumes a like.
         *  
         * @tparam object_type the type of the relevant object
         * @tparam candidate_type the type of the candidates in the list
         * @tparam links_type the type of the links the candidate is holding
         * 
         **/
        template <typename object_type, typename candidate_type, typename links_type>
        struct navigation_kernel
        {
            const object_type *on = nullptr;
            dvector<candidate_type> candidates = {};
            typename dvector<candidate_type>::iterator next = candidates.end();
            links_type links;
        };

        /** A navigation state object used to cache the information of the
         * current navigation stream.
         **/
        struct navigation_state
        {
            /// Kernel for the surfaces
            navigation_kernel<surface, intersection, surface_link> surface_kernel;
            /// Kernel the portals
            navigation_kernel<portal, intersection, portal_links> portal_kernel;

            // Volume navigation
            int volume_index = -1;
            /// Distance to next
            scalar distance_to_next = std::numeric_limits<scalar>::infinity();

            /// The detector in which we are moving
            detector_type detector;

            /// The inspector type of this navigation engine
            inspector_type inspector;
        };

        __plugin::cartesian2 cart2;
        __plugin::polar2 pol2;
        __plugin::cylindrical2 cyl2;

        /** Navigation status() call which established the current navigation information.
         * 
         * @param navigation is the navigation cache object
         * @param track is the track infromation
         * @param trust is a boolean that indicates if you can trust already the last estimate
         * 
         * @return a navigation status to be further decoded 
         **/
        navigation_status status(navigation_state &navigation, track<context> &track)
        {
          
            return e_unknown;
        }

        /** Target function of the navigator, finds the next candidates and set the distance to next
         * 
         * @tparam track is the track templated on the transform
         * 
         * @param navigation is the navigation cache
         * @param track is the current track information
         * 
         * @return a navigaiton status
         **/
        navigation_status target(navigation_state &navigation, track<context> &track)
        {

            return e_unknown;
        }
    };

} // namespace detray