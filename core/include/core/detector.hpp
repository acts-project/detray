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
#include "utils/indexing.hpp"
#include "utils/containers.hpp"
#include "utils/enumerate.hpp"
#include "tools/concentric_cylinder_intersector.hpp"

#include <exception>
#include <functional>
#include <cmath>
#include <climits>
#include <string>

#include <iostream>

namespace detray
{

    /** Detector representation in detray.
     * 
     * @tparam transform_type is the type of the underlying algebra transform
     */
    template <typename transform_type>
    class detector
    {

    public:
        /** Constructor with name 
         * @param name of the detector 
        */
        detector(const std::string &name = "unkown") : _name(name) {}

        // Algebra
        using point3 = typename transform_type::point3;
        using vector3 = typename transform_type::vector3;
        using point2 = __plugin::cartesian2::point2;
        using transform3 = transform_type;

        using surface_intersection = intersection<scalar, point3, point2>;

        // Indexing
        using typed_dindex = darray<dindex, 2>;
        using typed_dindex_range = dtuple<dindex, dindex_range>;

        // Surface finding function
        using local_object_finder = std::function<dvector<dindex>(const point2 &, const darray<unsigned int, 2> &)>;

        /** Nested volume class */
        struct volume
        {
            dvector<dindex> portal_surface_indices = {};
            dvector<dindex> surface_indices = {};
            dvector<dindex> surface_finder_indices = {};
            dindex index = 0;

            std::string name = "unknown";
            darray<scalar, 4> volume_bounds = {0,
                                               std::numeric_limits<scalar>::infinity(),
                                               -std::numeric_limits<scalar>::infinity(),
                                               std::numeric_limits<scalar>::infinity()};

            volume(const std::string &vname, const darray<scalar, 4> &bounds)
                : name(vname), volume_bounds(bounds)
            {
            }
        };

        /** Add a volume to the cdetector 
         * 
         * @param name the name of the new volume
         * @param bounds are the cylindrical bounds of the volume
         * 
         * @return the index of the new volume
         **/
        dindex new_volume(const std::string &name, const darray<scalar, 4> &bounds)
        {
            volume vol(name, bounds);
            vol.index = _volumes.size();
            _volumes.push_back(vol);
            return vol.index;
        }

        /** Portal section ******************************************************
         * 
         * Portals are masks that are applied to surfaces and then point to volumes
         * 
         * portal_links are [ opposite volume, along volume, opposite object finder, along volume finder ]
         * 
         **/
        using portal_links = darray<dindex, 4>;

        using portal_rectangle_mask = rectangle2<scalar, planar_intersector, portal_links>;
        using portal_rectangles = dvector<portal_rectangle_mask>;
        using protal_trapezoid_mask = trapezoid2<scalar, planar_intersector, portal_links>;
        using portal_trapezoids = dvector<protal_trapezoid_mask>;
        using portal_cylinder_mask = cylinder3<scalar, false, concentric_cylinder_intersector, portal_links>;
        using portal_cylinders = dvector<portal_cylinder_mask>;
        using portal_disc_mask = ring2<scalar, planar_intersector, portal_links>;
        using portal_discs = dvector<portal_disc_mask>;
        using portal_type_map = dmap<dindex, dindex>;
        using portal_surface = surface<dindex, typed_dindex_range, dindex>;

        /** Method to add a list of portals to the portal tuple
         * 
         * @tparam portal_type the type of the portal mask
         * 
         * @param volume_index the volume index to which this portal should be added
         * @param transform the transform of the portal surface
         * @param portals the vector of portal masks
         * 
         * @returns the index of the newly added portal surface
         **/
        template <typename portal_type>
        dindex
        add_portal_surface(dindex volume_index, transform_type &&transform, const dvector<portal_type> &portals)
        {
            auto &volume = _volumes[volume_index];
            // Get the boundary group, record the index and insert the portals
            auto &group = std::get<dvector<portal_type> >(_portal_masks);
            dindex index_start = group.size();
            dindex index_end = static_cast<dindex>(index_start + portals.size() - 1);
            group.insert(group.end(), portals.begin(), portals.end());
            // Create a range of portal masks
            dindex_range range = {index_start, index_end};
            dindex type = _portal_types.find(portal_type::mask_identifier)->second;
            // Record the portal index
            volume.portal_surface_indices.push_back(_portal_surfaces.size());
            typed_dindex_range links = {type, range};
            // Record the transform index
            dindex transform_index = _portal_transforms.size();
            _portal_transforms.push_back(std::move(transform));
            _portal_surfaces.push_back(portal_surface(std::move(transform_index), std::move(links), std::move(volume_index), false));
            return transform_index;
        }

        /** Update portal links
         *
         * @param portal_mask The portal mask which will be updated
         * @param additional_link The new link to be added to the portal
         *
         */
        template <typename mask_type>
        void update_portal_links(mask_type &mask, portal_links additional_link) noexcept(false)
        {
            auto &mask_link = mask.links();
            if (mask_link.size() != additional_link.size())
            {
                throw std::runtime_error("detray::update_protal_links(...) called with inconsistent link numbers.");
            }
            for (dindex il = 0; il < mask_link.size(); ++il)
            {
                auto &link = mask_link[il];
                auto &new_link = additional_link[il];
                if (link != dindex_invalid and new_link != dindex_invalid)
                {
                    std::string error_message = "detray::update_portal_links(...) is trying to overwrite valid link " + std::to_string(link) + std::string(" with ") + std::to_string(new_link);
                    throw std::runtime_error(error_message);
                }
                link = (link == dindex_invalid) ? new_link : link;
            }
        }

        /** Re-use a portal surface, but set a new indexing.
         * 
         * @param volume_index to which this portal belongs to (now as well)
         * @param portal_index The global portal index for this portal
         * @param updated_link The new link to be added to the portal
         * 
         */
        void reuse_portal_surface(dindex volume_index, dindex portal_index, portal_links additional_link)
        {
            auto &volume = _volumes[volume_index];
            auto &portal = _portal_surfaces[portal_index];

            volume.portal_surface_indices.push_back(portal_index);

            const auto &typed_mask_range = portal.mask();
            const auto &type_mask = std::get<0>(typed_mask_range);
            const auto &mask_range = std::get<1>(typed_mask_range);
            if (type_mask == 0)
            {
                auto &mask_group = std::get<0>(_portal_masks);
                for (dindex im = mask_range[0]; im <= mask_range[1]; ++im)
                {
                    update_portal_links(mask_group[im], additional_link);
                }
            }
            else if (type_mask == 1)
            {
                auto &mask_group = std::get<1>(_portal_masks);
                for (dindex im = mask_range[0]; im <= mask_range[1]; ++im)
                {
                    update_portal_links(mask_group[im], additional_link);
                }
            }
            else if (type_mask == 2)
            {
                auto &mask_group = std::get<2>(_portal_masks);
                for (dindex im = mask_range[0]; im <= mask_range[1]; ++im)
                {
                    update_portal_links(mask_group[im], additional_link);
                }
            }
            else if (type_mask == 3)
            {
                auto &mask_group = std::get<3>(_portal_masks);
                for (dindex im = mask_range[0]; im <= mask_range[1]; ++im)
                {
                    update_portal_links(mask_group[im], additional_link);
                }
            }
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
        using surface_links = bool;

        using surface_type_map = dmap<dindex, dindex>;
        using detector_surface = surface<dindex, typed_dindex_range, dindex>;

        /** Method to add a list of portals to the portal tuple
         * 
         * @tparam mask_type the type of the portal mask
         * 
         * @param volume_index the volume index to which these surfaces should be added
         * @param transform the transform of the portal surface
         * @param portals the vector of portal masks
         * 
         **/
        template <typename mask_type>
        void add_surfaces(dindex volume_index, dvector<transform_type> transforms, const typename mask_type::mask_values &mask_values)
        {
            auto &volume = _volumes[volume_index];
            // Get the boundary group, record the index and insert the portals
            auto &mask_group = std::get<dvector<mask_type> >(_surface_masks);
            dindex mask_index = mask_group.size();
            mask_type mask;
            mask = mask_values;
            mask_group.push_back(std::move(mask));
            dindex type = _surface_types.find(mask_type::mask_identifier)->second;
            typed_dindex_range typed_mask_range = {type, {mask_index, mask_index}};
            for (auto transform : transforms)
            {
                volume.surface_indices.push_back(_surfaces.size());
                dindex transform_index = _surface_transforms.size();
                _surface_transforms.push_back(std::move(transform));
                _surfaces.push_back(detector_surface(std::move(transform_index), std::move(typed_mask_range), std::move(volume_index), false));
            }
        }

        /** Add surface finders to a volume
         *
         * @param surface_finders the local finders that are to be added to this volume
         */
        dindex_range add_surface_finders(dvector<local_object_finder> surface_finders)
        {
            dindex finder_start_index = _surface_finders.size();
            _surface_finders.insert(_surface_finders.begin(), surface_finders.begin(), surface_finders.end());
            return {finder_start_index, finder_start_index + surface_finders.size()};
        }

        /** Const access method for the detector name */
        const std::string &name() const { return _name; }

        /** Const access for the volumes contained in this detector */
        const dvector<volume> &volumes() const { return _volumes; }

        /** Const access for all portal transforms  */
        const dvector<transform_type> &portal_transforms() const { return _portal_transforms; }

        /** Const access for all portal surfaces  */
        const dvector<portal_surface> &portal_surfaces() const { return _portal_surfaces; }

        /** Const access for all portal masks  */
        const auto &portal_masks() const { return _portal_masks; }

        /** Const access for all portal types */
        const portal_type_map portal_types() const { return _portal_types; }

        /** Const access for all surface transforms  */
        const dvector<transform_type> &surface_transforms() const { return _surface_transforms; }

        /** Const access for all surfaces */
        const dvector<detector_surface> &surfaces() const { return _surfaces; }

        /** Const access for all surface masks */
        const auto &surface_masks() const { return _surface_masks; }

        /** Const access for all surface mask types */
        const surface_type_map surface_types() const { return _surface_types; }

    private:
        std::string _name = "unknown";
        dvector<volume> _volumes;

        dvector<transform_type> _portal_transforms;
        dvector<portal_surface> _portal_surfaces;
        dtuple<portal_rectangles, portal_trapezoids, portal_cylinders, portal_discs> _portal_masks;
        portal_type_map _portal_types = {{rectangle_mask::mask_identifier, 0},
                                         {trapezoid_mask::mask_identifier, 1},
                                         {cylinder_mask::mask_identifier, 2},
                                         {disc_mask::mask_identifier, 3}};

        dvector<transform_type> _surface_transforms; //!< @todo change to contextual container
        dvector<detector_surface> _surfaces;
        dtuple<rectangles, trapezoids, cylinders, discs> _surface_masks;
        surface_type_map _surface_types = {{rectangle_mask::mask_identifier, 0},
                                           {trapezoid_mask::mask_identifier, 1},
                                           {cylinder_mask::mask_identifier, 2},
                                           {disc_mask::mask_identifier, 3}};

        dvector<local_object_finder> _surface_finders = {};
    };

} // namespace detray