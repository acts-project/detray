
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
namespace detray
{
    template <typename detector_type>
    struct navigator
    {
        using detector_surface = typename detector_type::detector_surface;
        using surface_intersection = typename detector_type::surface_intersection;
        using surface_links = typename detector_type::surface_links;
        using portal_surface = typename detector_type::portal_surface;
        using portal_links = typename detector_type::portal_links;
        using transform_type = typename detector_type::transform3;

        // Indexing
        using guaranteed_index = typename detector_type::guaranteed_index;

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
            // For the surfaces
            navigation_kernel<detector_surface, surface_intersection, surface_links> surfaces;
            // For the portals
            navigation_kernel<portal_surface, surface_intersection, portal_links> portals;
            // Volume navigation
            int volume_index = -1;
            // Distance to next
            scalar distance_to_next = std::numeric_limits<scalar>::infinity();

            // The detector in which we are moving
            detector_type detector;
        };

        __plugin::cartesian2 cart2;
        __plugin::polar2 pol2;
        __plugin::cylindrical2 cyl2;

        /** Specialized method to update an intersection when the mask group is resolved
         * 
         * @tparam links_type is the type of links object held for this mask group
         * @tparam surface_type is the type of the transfrom used
         * @tparam local_type  is the type of the local cast operation
         * @tparam mask_group is the type of the split out mask group from all masks
         * @tparam mask_range is the type of the according mask range object
         * 
         * @param sfi is the surface intersection to be updated
         * @param links is the associated object link to be updated
         * @param track is the track information 
         * @param surface is the surface to be intersected
         * @param local is the local coordiante cast
         * @param masks is the associated (and split out) mask group
         * @param range is the range list of masks to be processed
         * 
         * @note since all masks are on the same surface, only maximally one mask solution at a time
         * is possbile.
         * 
         * @return a boolean indicating if the update succeeded.
         **/
        template <typename links_type, typename local_type, typename mask_group, typename mask_range>
        bool update_intersection_by_mask(surface_intersection &sfi,
                                         links_type &links,
                                         track<transform_type> &track,
                                         const transform_type &trf,
                                         const local_type &local,
                                         const mask_group &masks,
                                         const mask_range &range)
        {
            for (guaranteed_index i = range[0]; i <= range[1]; ++i)
            {
                auto &mask = masks[i];
                auto is = mask.intersector();
                sfi = is.intersect(trf, track, local, mask);
                if (sfi.status == e_inside)
                {
                    links = mask.links();
                    return true;
                }
            }
            return false;
        }

        /** Actual kernel method that updates the updates the intersections
         * 
         * @tparam surface_type the type of surfaces to be processed
         * @tparam links_type the type of the links held by the masks
         * @tparam mask_type_map a mapping between indices and mask columsn in the tuple
         * @tparam mask_container the type of the tupled mask container
         * 
         * @param sfi is the surface intersection to be updated
         * @param links is the associated object link to be updated
         * @param track is the track information 
         * @param trf is the transform of surface 
         * @param surface  is the surface to be intersected
         * @param mask_types is the mapping between mask types and tuple entry
         * @param masks is the associated (and split out) mask group
         * @param trust is a boolean flag whether you can trust the last estimate already
         *
         * @return a boolean indicating if the update was successful
         **/
        template <typename surface_type, typename links_type, typename mask_type_map, typename mask_container>
        bool update_intersection(surface_intersection &sfi,
                                 links_type &links,
                                 track<transform_type> &track,
                                 const transform_type &trf,
                                 const surface_type &surface,
                                 const mask_type_map &mask_types,
                                 const mask_container &masks,
                                 bool trust = false)
        {
            // That the guaranteed index in the container
            const auto &typed_mask_range = surface.mask();
            if (std::get<0>(typed_mask_range) == 0)
            {
                const auto &mask_group = std::get<0>(masks);
                return update_intersection_by_mask(sfi, links, track, trf, cart2, mask_group, std::get<1>(typed_mask_range));
            }
            else if (std::get<0>(typed_mask_range) == 1)
            {
                const auto &mask_group = std::get<1>(masks);
                return update_intersection_by_mask(sfi, links, track, trf, cart2, mask_group, std::get<1>(typed_mask_range));
            }
            else if (std::get<0>(typed_mask_range) == 2)
            {
                const auto &mask_group = std::get<2>(masks);
                return update_intersection_by_mask(sfi, links, track, trf, cyl2, mask_group, std::get<1>(typed_mask_range));
            }
            else if (std::get<0>(typed_mask_range) == 3)
            {
                const auto &mask_group = std::get<3>(masks);
                return update_intersection_by_mask(sfi, links, track, trf, cart2, mask_group, std::get<1>(typed_mask_range));
            }
            return false;
        }

        /** Navigation status() call which established the current navigation information.
         * 
         * @param navigation is the navigation cache object
         * @param track is the track infromation
         * @param trust is a boolean that indicates if you can trust already the last estimate
         * 
         * @return a navigation status to be further decoded 
         **/
        navigation_status status(navigation_state &navigation, track<transform_type> &track, bool trust = false)
        {
            if (trust)
            {
                const auto &surfaces = navigation.detector.surfaces();
                const auto &portal_surfaces = navigation.detector.portal_surfaces();

                // Trust the stepper towards the internal surface
                if (navigation.surfaces.next != navigation.surfaces.candidates.end())
                {
                    navigation.surfaces.on = &(surfaces[(*navigation.surfaces.next).index]);
                    ++navigation.surfaces.next;
                    return e_on_surface;
                }

                // Trust the stepper towards the portal surface
                if (navigation.portals.next != navigation.portals.candidates.end())
                {
                    navigation.portals.on = &(portal_surfaces[(*navigation.portals.next).index]);
                    auto attacheddirection = (*navigation.portals.next).direction;
                    navigation.volume_index = navigation.portals.links[attacheddirection];
                    navigation.portals.candidates.clear();
                    return e_on_portal;
                }
            }
            return e_unknown;
        }

        /** Target function of the navigator, finds the next candidates and set the distance to next
         * 
         * @tparam track is the track templated on the transform
         * 
         * @param navigation is the navigation cache
         * @param track is the current track information
         * @param trust is a flag indicating if you can trust the already found intersections
         * 
         * @return a navigaiton status
         **/
        navigation_status target(navigation_state &navigation, track<transform_type> &track, bool trust = false)
        {

            const auto &volumes = navigation.detector.volumes();

            if (navigation.volume_index >= 0)
            {
                // Volume information
                const auto &volume = volumes[navigation.volume_index];
                // Fresh volume call
                if (navigation.surfaces.candidates.empty() and navigation.portals.candidates.empty())
                {
                    const auto &surfaces = navigation.detector.surfaces();
                    const auto &surface_transforms = navigation.detector.surface_transforms();
                    const auto &surface_types = navigation.detector.surface_types();
                    const auto &surface_masts = navigation.detector.surface_masks();

                    // This is the code without surface_finder (start version for the moment)
                    navigation.surfaces.candidates.reserve(volume.surface_indices.size());
                    for (auto isi : volume.surface_indices)
                    {
                        surface_intersection sfi;
                        sfi.index = isi;
                        const auto &surface = surfaces[isi];
                        // @todo this will become a contextual call
                        const auto &surface_transform = surface_transforms[surface.transform()];
                        update_intersection(sfi, navigation.surfaces.links, track, surface_transform, surface, surface_types, surface_masts);
                        if (sfi.status == e_inside and sfi.path > track.overstep_tolerance)
                        {
                            navigation.surfaces.candidates.push_back(std::move(sfi));
                        }
                    }
                    if (not navigation.surfaces.candidates.empty())
                    {
                        std::sort(navigation.surfaces.candidates.begin(), navigation.surfaces.candidates.end());
                        navigation.surfaces.next = navigation.surfaces.candidates.begin();
                        navigation.distance_to_next = (*navigation.surfaces.next).path;
                        // we need to sort them still
                        return e_towards_surface;
                    }
                }
                else if (navigation.surfaces.next == navigation.surfaces.candidates.end())
                {
                    navigation.surfaces.candidates.clear();
                    navigation.surfaces.next = navigation.surfaces.candidates.end();

                    const auto &portal_transforms = navigation.detector.portal_transforms();
                    const auto &portal_surfaces = navigation.detector.portal_surfaces();
                    const auto &portal_types = navigation.detector.portal_types();
                    const auto &portal_masks = navigation.detector.portal_masks();

                    // This is the code without portal_finder (start verison for the moment)
                    navigation.portals.candidates.reserve(volume.portal_surface_indices.size());
                    bool portal_hit = false;
                    for (auto psi : volume.portal_surface_indices)
                    {
                        surface_intersection sfi;
                        sfi.index = psi;
                        if (not portal_hit)
                        {
                            const auto &portal_surface = portal_surfaces[sfi.index];
                            const auto &portal_transform = portal_transforms[portal_surface.transform()];
                            update_intersection(sfi, navigation.portals.links, track, portal_transform, portal_surface, portal_types, portal_masks);
                            portal_hit = (sfi.status == e_inside and sfi.path > track.overstep_tolerance);
                            if (portal_hit)
                            {
                                // That's not the full story, we should keep all boundaries
                                navigation.portals.candidates.push_back(std::move(sfi));
                            }
                        }
                    }
                    navigation.portals.next = navigation.portals.candidates.begin();
                    navigation.distance_to_next = (*navigation.portals.next).path;
                    return e_towards_portal;
                }
            }
            return e_unknown;
        }
    };

} // namespace detray