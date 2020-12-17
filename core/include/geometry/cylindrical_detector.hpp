/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "core/surface.hpp"
#include "masks/rectangle2.hpp"
#include "masks/trapezoid2.hpp"
#include "masks/cylinder3.hpp"
#include "masks/ring2.hpp"
#include "utils/containers.hpp"
#include "tools/concentric_cylinder_intersector.hpp"

#include <cmath>
#include <climits>
#include <iostream>

namespace detray
{

    template <typename transform_type>
    class cylindrical_detector
    {

    public:
        // Algebra
        using point3 = typename transform_type::point3;
        using vector3 = typename transform_type::vector3;
        using context = typename transform_type::context;
        using point2 = __plugin::cartesian2::point2;

        using surface_intersection = intersection<scalar, point3, point2>;

        // Indexing
        using optional_index = int;
        using guaranteed_index = unsigned long;
        using typed_guaranteed_index = darray<guaranteed_index, 2>;
        using guaranteed_range = darray<guaranteed_index, 2>;
        using typed_guaranteed_range = dtuple<guaranteed_index, guaranteed_range>;

        // Surface finding function
        using surface_finder = std::function<const guaranteed_range &(const point3 &, const vector3 &, const vector3 &, scalar)>;

        /** Nested volume class */
        class volume
        {

            /*** The cylindrical detector is declared friend to call the constructor */
            friend class cylindrical_detector<transform_type>;

        private:
            dvector<guaranteed_index> _portal_surface_indices = {};
            dvector<surface_finder> _portal_surface_finder = {};
            dvector<guaranteed_index> _internal_surface_indices = {};
            dvector<surface_finder> _internal_surface_finder = {};
            guaranteed_index _volume_index = 0;

            darray<scalar, 4> _volume_bounds = {0,
                                                std::numeric_limits<scalar>::infinity(),
                                                -std::numeric_limits<scalar>::infinity(),
                                                std::numeric_limits<scalar>::infinity()};

            volume(const darray<scalar, 4> &bounds)
                : _volume_bounds(bounds)
            {
            }
        };

        /** Volume handling section ******************************************************
         * 
         * Volumes are a container of boundary surface indices and contained surface indices
         **/
    private:
        dvector<volume> _volumes;

    public:
        /** Add a volume to the detector 
         * 
         * @tparam bounds are the cylindrical bounds of the volume
         **/
        volume &new_volume(const darray<scalar, 4> &bounds)
        {
            volume vol(bounds);
            vol._volume_index = _volumes.size();
            _volumes.push_back(vol);
            return _volumes[vol._volume_index];
        }

        /** Portal section ******************************************************
         * 
         * Portals are masks that are applied to surfaces and then point to volumes
         * 
         **/
        using portal_destinations = darray<optional_index, 2>;

        using portal_rectangle_mask = rectangle2<scalar, planar_intersector, portal_destinations>;
        using portal_rectangles = dvector<portal_rectangle_mask>;
        using protal_trapezoid_mask = trapezoid2<scalar, planar_intersector, portal_destinations>;
        using portal_trapezoids = dvector<protal_trapezoid_mask>;
        using portal_cylinder_mask = cylinder3<scalar, false, concentric_cylinder_intersector, portal_destinations>;
        using portal_cylinders = dvector<portal_cylinder_mask>;
        using portal_disc_mask = ring2<scalar, planar_intersector, portal_destinations>;
        using portal_discs = dvector<portal_disc_mask>;
        using portal_type_map = dmap<guaranteed_index, guaranteed_index>;
        using portal_surface = surface<transform_type, typed_guaranteed_range, guaranteed_index>;

        // Member & method section: for portals
    private:
        dvector<portal_surface> _portal_surfaces;
        dtuple<portal_rectangles, portal_trapezoids, portal_cylinders, portal_discs> _portals;
        portal_type_map _portal_types = {{rectangle_mask::mask_identifier, 0},
                                         {trapezoid_mask::mask_identifier, 1},
                                         {cylinder_mask::mask_identifier, 2},
                                         {disc_mask::mask_identifier, 3}};

    public:
        /** Method to add a list of portals to the portal tuple
         * 
         * @tparam portal_type the type of the portal mask
         * 
         * @param transform the transform of the portal surface
         * @param portals the vector of portal masks
         * @param volume to which this portal belongs to
         **/
        template <typename portal_type>
        void add_portal_surface(transform_type &&transform, const dvector<portal_type> &portals, volume &volume)
        {
            // Get the boundary group, record the index and insert the portals
            auto &group = std::get<dvector<portal_type>>(_portals);
            guaranteed_index index_start = group.size();
            guaranteed_index index_end = static_cast<guaranteed_index>(index_start + portals.size() - 1);
            group.insert(group.end(), portals.begin(), portals.end());
            // Create a range of portal masks
            guaranteed_range range = {index_start, index_end};
            guaranteed_index type = _portal_types.find(portal_type::mask_identifier)->second;
            // Record the portal index
            volume._portal_surface_indices.push_back(_portal_surfaces.size());
            guaranteed_index volume_index = volume._volume_index;
            typed_guaranteed_range links = {type, range};
            _portal_surfaces.push_back(portal_surface(std::move(transform), std::move(links), std::move(volume_index), false));
        }

        /** Internal surface section ***********************************************
         * 
         * Internal surfaces are all kind of navigation surfaces within a volume
         * 
         **/
        using rectangle_mask = rectangle2<scalar>;
        using trapezoid_mask = trapezoid2<scalar>;
        using cylinder_mask = cylinder3<scalar, false, concentric_cylinder_intersector>;
        using disc_mask = ring2<scalar, planar_intersector>;
        using rectangles = dvector<rectangle_mask>;
        using trapezoids = dvector<trapezoid_mask>;
        using cylinders = dvector<cylinder_mask>;
        using discs = dvector<disc_mask>;

        using surface_type_map = dmap<guaranteed_index, guaranteed_index>;
        using internal_surface = surface<transform_type, typed_guaranteed_range, guaranteed_index>;

        // Member & method section: for portals
    private:
        dvector<internal_surface> _internal_surfaces;
        dtuple<rectangles, trapezoids, cylinders, discs> _internal_masks;
        surface_type_map _internal_types = {{rectangle_mask::mask_identifier, 0},
                                            {trapezoid_mask::mask_identifier, 1},
                                            {cylinder_mask::mask_identifier, 2},
                                            {disc_mask::mask_identifier, 3}};

    public:
        /** Method to add a list of portals to the portal tuple
         * 
         * @tparam portal_type the type of the portal mask
         * 
         * @param transform the transform of the portal surface
         * @param portals the vector of portal masks
         * @param volume to which this portal belongs to
         **/
        template <typename mask_type>
        void add_internal_surfaces(dvector<transform_type> transforms, const typename mask_type::mask_values &mask_values, volume &volume)
        {
            // Get the boundary group, record the index and insert the portals
            auto &mask_group = std::get<dvector<mask_type>>(_internal_masks);
            guaranteed_index mask_index = mask_group.size();
            mask_type mask;
            mask = mask_values;
            mask_group.push_back(std::move(mask));
            guaranteed_index type = _internal_types.find(mask_type::mask_identifier)->second;
            typed_guaranteed_range typed_mask_range = {type, {mask_index, mask_index}};
            // Record the portal index
            volume._internal_surface_indices.push_back(_internal_surfaces.size());
            guaranteed_index volume_index = volume._volume_index;
            for (auto transform : transforms)
            {
                _internal_surfaces.push_back(internal_surface(std::move(transform), std::move(typed_mask_range), std::move(volume_index), false));
            }
        }

        // ----------- Navigation section --------
        enum navigation_status : int
        {
            e_unknown = -1,
            e_towards_surface = 0,
            e_on_surface = 1,
            e_towards_portal = 2,
            e_on_portal = 3,
        };

        struct navigation_state
        {
            // Internal surface navigation
            const internal_surface *surface = nullptr;
            dvector<surface_intersection> surface_candidates = {};
            typename dvector<surface_intersection>::iterator surface_iterator = surface_candidates.end();
            // Portal navigation
            const portal_surface *portal = nullptr;
            dvector<surface_intersection> portal_candidates = {};
            typename dvector<surface_intersection>::iterator portal_iterator = portal_candidates.end();
            // Distance to next
            scalar distance_to_next = std::numeric_limits<scalar>::infinity();
            // Volume navigation
            optional_index volume_index = -1;
        };

        __plugin::cartesian2 cart2;
        __plugin::polar2 pol2;
        __plugin::cylindrical2 cyl2;

        template <typename surface_type, typename local_type, typename mask_group, typename mask_range>
        bool update_intersection_by_mask(surface_intersection &sfi, track<transform_type> &track, const surface_type &surface, const local_type &local, const mask_group &masks, const mask_range &range)
        {
            // std::cout << "--- retrieved mask range is " << range[0] << ", " << range[1] << std::endl;
            for (guaranteed_index i = range[0]; i <= range[1]; ++i)
            {
                // std::cout << "---- attempt for mask " << i << std::endl;
                auto &mask = masks[i];
                auto is = mask.intersector();
                sfi = is.intersect(surface, track, local, mask);
                if (sfi._status == e_inside)
                {
                    return true;
                }
            }
            return false;
        }

        template <typename surface_container, typename mask_type_map, typename mask_container>
        bool update_intersection(surface_intersection &sfi, track<transform_type> &track,
                                 const surface_container &surfaces,
                                 const mask_type_map &mask_types,
                                 const mask_container &masks,
                                 bool trust = false)
        {
            // That the guaranteed index in the container
            const auto &surface = surfaces[sfi._index];
            const auto &typed_mask_range = surface.mask();
            if (std::get<0>(typed_mask_range) == 0)
            {
                // std::cout << "-- trying rectangle -- " << std::endl;
                const auto &mask_group = std::get<0>(masks);
                return update_intersection_by_mask(sfi, track, surface, cart2, mask_group, std::get<1>(typed_mask_range));
            }
            else if (std::get<0>(typed_mask_range) == 1)
            {
                //std::cout << "-- trying trapezoid -- " << std::endl;
                const auto &mask_group = std::get<1>(masks);
                return update_intersection_by_mask(sfi, track, surface, cart2, mask_group, std::get<1>(typed_mask_range));
            }
            else if (std::get<0>(typed_mask_range) == 2)
            {
                //std::cout << "-- trying cylinder -- " << std::endl;
                const auto &mask_group = std::get<2>(masks);
                return update_intersection_by_mask(sfi, track, surface, cyl2, mask_group, std::get<1>(typed_mask_range));
            }
            else if (std::get<0>(typed_mask_range) == 3)
            {
                //std::cout << "-- trying disc -- " << std::endl;
                const auto &mask_group = std::get<3>(masks);
                return update_intersection_by_mask(sfi, track, surface, cart2, mask_group, std::get<1>(typed_mask_range));
            }
            return false;
        }

        navigation_status status(navigation_state &navigation, track<transform_type> &track, bool trust = false)
        {
            //std::cout << "--> status() " << std::endl;
            if (trust)
            {
                // Trust the stepper towards the internal surface
                if (navigation.surface_iterator != navigation.surface_candidates.end())
                {
                    navigation.surface = &(_internal_surfaces[(*navigation.surface_iterator)._index]);
                    // std::cout << "--- On surface, trying next ... " << std::endl;
                    ++navigation.surface_iterator;
                    return e_on_surface;
                }

                // Trust the stepper towards the portal surface
                if (navigation.portal_iterator != navigation.portal_candidates.end())
                {
                    navigation.portal = &(_portal_surfaces[(*navigation.portal_iterator)._index]);   
                    // std::cout << "--- On portal, moving to next volume " << std::endl;                 
                    navigation.portal_candidates.clear();
                    return e_on_portal;
                }
            }
            return e_unknown;
        }

        // Navigation target function
        navigation_status target(navigation_state &navigation, track<transform_type> &track, bool trust = false)
        {
            // std::cout << "--> target() " << std::endl;
            if (navigation.volume_index >= 0)
            {
                // Volume information
                const auto &volume = _volumes[navigation.volume_index];
                // Fresh volume call
                if (navigation.surface_candidates.empty() and navigation.portal_candidates.empty())
                {
                    // std::cout << "--- Fresh volume call --- " << std::endl;
                    // This is the code without surface_finder (start version for the moment)
                    navigation.surface_candidates.reserve(volume._internal_surface_indices.size());
                    for (auto isi : volume._internal_surface_indices)
                    {
                        surface_intersection sfi;
                        sfi._index = isi;
                        update_intersection(sfi, track, _internal_surfaces, _internal_types, _internal_masks);
                        if (sfi._status == e_inside and sfi._path > track.overstep_tolerance)
                        {
                            navigation.surface_candidates.push_back(std::move(sfi));
                        }
                    }
                    // std::cout << "--- Found " << navigation.surface_candidates.size() << " surface candidates " << std::endl;
                    if (not navigation.surface_candidates.empty())
                    {
                        std::sort(navigation.surface_candidates.begin(), navigation.surface_candidates.end());
                        navigation.surface_iterator = navigation.surface_candidates.begin();
                        navigation.distance_to_next = (*navigation.surface_iterator)._path;
                        // we need to sort them still
                        // std::cout << "--- targetting surface candidate  " << std::endl;
                        return e_towards_surface;
                    }
                }
                else if (navigation.surface_iterator == navigation.surface_candidates.end())
                {
                    navigation.surface_candidates.clear();
                    navigation.surface_iterator = navigation.surface_candidates.end();
                    // std::cout << "--- No more surface candidates, get portal candidates" << std::endl;
                    // This is the code without portal_finder (start verison for the moment)
                    navigation.portal_candidates.reserve(volume._portal_surface_indices.size());
                    bool portal_hit = false;
                    for (auto psi : volume._portal_surface_indices)
                    {
                        surface_intersection sfi;
                        sfi._index = psi;
                        if (not portal_hit)
                        {
                            update_intersection(sfi, track, _portal_surfaces, _portal_types, _portals);
                            portal_hit = (sfi._status == e_inside and sfi._path > track.overstep_tolerance);
                            if (portal_hit){
                                // That's not the full story, we should keep all boundaries
                                navigation.portal_candidates.push_back(std::move(sfi));
                            }
                        }
                    }
                    navigation.portal_iterator = navigation.portal_candidates.begin();
                    navigation.distance_to_next = (*navigation.portal_iterator)._path;
                    // std::cout << "--- targetting portal candidate " << std::endl;
                    return e_towards_portal;
                }
            }
            return e_unknown;
        }
    };

} // namespace detray