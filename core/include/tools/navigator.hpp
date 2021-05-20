
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/intersection.hpp"
#include "core/track.hpp"
#include "utils/indexing.hpp"
#include "utils/enumerate.hpp"
#include "tools/intersection_kernel.hpp"

#include <iostream>

namespace detray
{

    /** A void inpector that does nothing.
     * 
     * Inspectors can be plugged in to understand the
     * the current navigation state information.
     * 
     */
    struct void_inspector
    {
        template <typename state_type>
        void operator()(const state_type & /*ignored*/) {}
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

        /** Navigation status flag */
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

        /** Navigation trust level */
        enum navigation_trust_level : int
        {
            e_not = 0,  // re-evalute the candidates all over
            e_fair = 1, // re-evaluate the distance & order of the (preselected) candidates
            e_high = 3, // re-evaluate the distance to the next candidate
            e_full = 4  // trust fully: distance & next candidate
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
            dindex volume_index = dindex_invalid;
            /// Distance to next
            scalar distance_to_next = std::numeric_limits<scalar>::infinity();
            /// The on surface tolerance
            scalar on_surface_tolerance = 1e-5;

            /// The inspector type of this navigation engine
            inspector_type inspector;

            /// The navigation status
            navigation_status status = e_unknown;
            /// The navigation trust level
            navigation_trust_level trust_level = e_not;
        };

        __plugin::cartesian2 cart2;
        __plugin::polar2 pol2;
        __plugin::cylindrical2 cyl2;

        /// The detector in which we are moving
        detector_type detector;

        /** Constructor with move constructor
         * 
         * @param d the detector for this navigator
         */
        navigator(detector_type &&d) : detector(std::move(d)) {}

        /** Navigation status() call which established the current navigation information.
         * 
         * @param navigation is the navigation cache object
         * @param track is the track infromation
         * @param trust is a boolean that indicates if you can trust already the last estimate
         * 
         * @return a navigation status to be further decoded 
         **/
        void status(navigation_state &navigation, const track<context> &track)
        {
            // Retrieve the volume
            const auto &volume = detector.indexed_volume(track.pos);
            navigation.volume_index = volume.index();

            // The navigation is not yet initialized, doing it now
            if (navigation.volume_index == dindex_invalid or navigation.trust_level == e_not)
            {
                // First try to get the surface candidates
                if (not volume.empty())
                {
                    initialize_kernel(navigation.surface_kernel, track, volume.surfaces());
                    sort_and_set(navigation, navigation.surface_kernel);
                }
                // If no surfaces are to processed, initialize the portals
                if (navigation.surface_kernel.candidates.empty())
                {
                    initialize_kernel(navigation.portal_kernel, track, volume.portals());
                    sort_and_set(navigation, navigation.portal_kernel);
                }
                return;
            }

            if (navigation.trust_level == e_high)
            {
                // Block to update current candidate
                if (navigation.surface_kernel.next != navigation.surface_kernel.candidates.end())
                {
                    // Only update the last intersection
                    dindex si = navigation.surface_kernel.next->index;
                    const auto &tfs = volume.surfaces().transforms();
                    const auto &msks = volume.surfaces().masks();
                    const auto &s = volume.surfaces().indexed_object(si);
                    auto [sfi, link] = intersect(track, s, tfs, msks);
                    sfi.index = si;
                    if (sfi.status == e_inside)
                    {
                        // Update the intersection
                        (*navigation.surface_kernel.next) = sfi;
                        navigation.distance_to_next = sfi.path;
                        // Trust fully again
                        navigation.trust_level = e_full;
                        return;
                    }
                }
            }

            return;
        }

        /** Target function of the navigator, finds the next candidates and set the distance to next
         * 
         * @param navigation is the navigation cache
         * @param track is the current track information
         * 
         * @return a navigaiton status
         **/
        void target(navigation_state &navigation, const track<context> &track)
        {
            // Retrieve the volume
            const auto &volume = detector.indexed_volume(track.pos);
            navigation.volume_index = volume.index();

            // Full trust level from that status call
            if (navigation.trust_level == e_full)
            {
                return;
            }
            return;
        }

        /** Helper method to intersect all objects of a constituents store
         * 
         * @tparam kernel_t the type of the kernel: surfaces/portals
         * @tparam constituents the type of the associated constituents
         * 
         * @param kernel [in,out] the kernel to be checked/initialized 
         * @param track the track information 
         * @param constituens the constituents to be checked
         * 
         */
        template <typename kernel_t, typename constituents_t>
        void initialize_kernel(kernel_t &kernel,
                               const track<context> &track,
                               const constituents_t &constituents)
        {

            // Get the number of candidates & run them throuth the kernel
            size_t n_surfaces = constituents.objects().size();
            kernel.candidates.reserve(n_surfaces);
            const auto &transforms = constituents.transforms();
            const auto &masks = constituents.masks();
            // Loop over all indexed surfaces, intersect and fill
            for (auto si : sequence({0, n_surfaces - 1}))
            {
                const auto &object = constituents.indexed_object(si);
                auto [sfi, link] = intersect(track, object, transforms, masks);
                sfi.index = si;
                if (sfi.status == e_inside)
                {
                    kernel.candidates.push_back(sfi);
                }
            }
        }

        /** Helper method to check and set
         *
         * @param navigation is the state for setting the distance to next
         * @param kernel is the kernel which should be updated
         */
        template <typename kernel_t>
        void check_and_set(navigation_state &navigation, kernel_t &kernel)
        {
            // Get the type of the kernel via a const expression
            constexpr bool kType = (std::is_same_v<kernel_t, navigation_kernel<surface, intersection, surface_link>>);
        }

        /** Helper method to sort within the kernel 
         *
         * @param navigation is the state for setting the distance to next
         * @param kernel is the kernel which should be updated
         */
        template <typename kernel_t>
        void sort_and_set(navigation_state &navigation, kernel_t &kernel)
        {

            // Get the type of the kernel via a const expression
            constexpr bool kType = (std::is_same_v<kernel_t, navigation_kernel<surface, intersection, surface_link>>);

            // Sort and set distance to next
            if (not kernel.candidates.empty())
            {
                navigation.trust_level = e_full;
                // Sort and set
                std::sort(kernel.candidates.begin(), kernel.candidates.end());
                kernel.next = kernel.candidates.begin();
                navigation.distance_to_next = kernel.next->path;
                if (navigation.distance_to_next < navigation.on_surface_tolerance)
                {
                    navigation.status = kType ? e_on_surface : e_on_portal;
                }
                navigation.status = kType ? e_towards_surface : e_towards_portal;
            }
        }
    };

} // namespace detray