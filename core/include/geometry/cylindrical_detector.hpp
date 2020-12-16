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

namespace detray
{

    template <typename transform_type>
    class cylindrical_detector
    {

    public:
        // Algebra
        using point3 = typename transform_type::point3;
        using vector3 = typename transform_type::vector3;

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
        using portal_cylinder_mask = cylinder3<scalar, false, concentric_cylinder_intersector, portal_destinations>;
        using portal_cylinders = dvector<portal_cylinder_mask>;
        using portal_disc_mask = ring2<scalar, planar_intersector, portal_destinations>;
        using portal_discs = dvector<portal_disc_mask>;
        using portal_type_map = dmap<guaranteed_index, guaranteed_index>;
        using portal_surface = surface<transform_type, typed_guaranteed_range, guaranteed_index>;

        // Member & method section: for portals
    private:
        dvector<portal_surface> _portal_surfaces;
        dtuple<portal_cylinders, portal_discs> _portals;
        portal_type_map _portal_types = {{portal_cylinder_mask::mask_identifier, 0}, {portal_disc_mask::mask_identifier, 1}};

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
            group.insert(group.end(), portals.begin(), portals.end());
            // Create a range of portal masks
            guaranteed_range range = {index_start, index_start + portals.size()};
            guaranteed_index type = _portal_types.find(portal_type::mask_identifier)->first;
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
        using inernal_surface = surface<transform_type, typed_guaranteed_index, guaranteed_index>;

        // Member & method section: for portals
    private:
        dvector<inernal_surface> _internal_surfaces;
        dtuple<rectangles, trapezoids, cylinders, discs> _internal_masks;
        surface_type_map _internal_types = {{rectangle_mask::mask_identifier, 0}, {trapezoid_mask::mask_identifier, 1}, {cylinder_mask::mask_identifier, 0}, {disc_mask::mask_identifier, 1}};

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
            guaranteed_index type = _portal_types.find(mask_type::mask_identifier)->first;
            typed_guaranteed_index typed_mask_index = {type, mask_index};
            // Record the portal index
            volume._internal_surface_indices.push_back(_internal_surfaces.size());
            guaranteed_index volume_index = volume._volume_index;
            for (auto transform : transforms)
            {
                _internal_surfaces.push_back(inernal_surface(std::move(transform), std::move(typed_mask_index), std::move(volume_index), false));
            }
        }
    };

} // namespace detray