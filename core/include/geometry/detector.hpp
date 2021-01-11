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
#include "tools/concentric_cylinder_intersector.hpp"

#include <functional>
#include <cmath>
#include <climits>
#include <string>

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

        /** Constructor with name */
        detector(const std::string& name = "unkown") :
         _name(name){}

        // Algebra
        using point3 = typename transform_type::point3;
        using vector3 = typename transform_type::vector3;
        using point2 = __plugin::cartesian2::point2;
        using transform3 = transform_type;

        using surface_intersection = intersection<scalar, point3, point2>;

        // Indexing
        using typed_guaranteed_index = darray<guaranteed_index, 2>;
        using typed_guaranteed_range = dtuple<guaranteed_index, guaranteed_range>;

        // Surface finding function
        using object_finder = std::function<const dvector<guaranteed_index> &(const point3 &, const vector3 &, const vector3 &, scalar)>;

        /** Nested volume class */
        struct volume
        {
            dvector<guaranteed_index> portal_surface_indices = {};
            dvector<object_finder> portal_surface_finder = {};
            dvector<guaranteed_index> surface_indices = {};
            dvector<object_finder> surface_finder = {};
            guaranteed_index index = 0;

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
         * @tparam bounds are the cylindrical bounds of the volume
         **/
        volume &new_volume(const std::string &name, const darray<scalar, 4> &bounds)
        {
            volume vol(name, bounds);
            vol.index = _volumes.size();
            _volumes.push_back(vol);
            return _volumes[vol.index];
        }

        /** Portal section ******************************************************
         * 
         * Portals are masks that are applied to surfaces and then point to volumes
         * 
         **/
        using portal_links = darray<optional_index, 2>;

        using portal_rectangle_mask = rectangle2<scalar, planar_intersector, portal_links>;
        using portal_rectangles = dvector<portal_rectangle_mask>;
        using protal_trapezoid_mask = trapezoid2<scalar, planar_intersector, portal_links>;
        using portal_trapezoids = dvector<protal_trapezoid_mask>;
        using portal_cylinder_mask = cylinder3<scalar, false, concentric_cylinder_intersector, portal_links>;
        using portal_cylinders = dvector<portal_cylinder_mask>;
        using portal_disc_mask = ring2<scalar, planar_intersector, portal_links>;
        using portal_discs = dvector<portal_disc_mask>;
        using portal_type_map = dmap<guaranteed_index, guaranteed_index>;
        using portal_surface = surface<guaranteed_index, typed_guaranteed_range, guaranteed_index>;

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
            auto &group = std::get<dvector<portal_type>>(_portal_masks);
            guaranteed_index index_start = group.size();
            guaranteed_index index_end = static_cast<guaranteed_index>(index_start + portals.size() - 1);
            group.insert(group.end(), portals.begin(), portals.end());
            // Create a range of portal masks
            guaranteed_range range = {index_start, index_end};
            guaranteed_index type = _portal_types.find(portal_type::mask_identifier)->second;
            // Record the portal index
            volume.portal_surface_indices.push_back(_portal_surfaces.size());
            guaranteed_index volume_index = volume.index;
            typed_guaranteed_range links = {type, range};
            // Record the transform index
            guaranteed_index transform_index = _portal_transforms.size();
            _portal_transforms.push_back(std::move(transform));
            _portal_surfaces.push_back(portal_surface(std::move(transform_index), std::move(links), std::move(volume_index), false));
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

        using surface_type_map = dmap<guaranteed_index, guaranteed_index>;
        using detector_surface = surface<guaranteed_index, typed_guaranteed_range, guaranteed_index>;

        /** Method to add a list of portals to the portal tuple
         * 
         * @tparam portal_type the type of the portal mask
         * 
         * @param transform the transform of the portal surface
         * @param portals the vector of portal masks
         * @param volume to which this portal belongs to
         **/
        template <typename mask_type>
        void add_surfaces(dvector<transform_type> transforms, const typename mask_type::mask_values &mask_values, volume &volume)
        {
            // Get the boundary group, record the index and insert the portals
            auto &mask_group = std::get<dvector<mask_type>>(_surface_masks);
            guaranteed_index mask_index = mask_group.size();
            mask_type mask;
            mask = mask_values;
            mask_group.push_back(std::move(mask));
            guaranteed_index type = _surface_types.find(mask_type::mask_identifier)->second;
            typed_guaranteed_range typed_mask_range = {type, {mask_index, mask_index}};
            guaranteed_index volume_index = volume.index;
            for (auto transform : transforms)
            {
                volume.surface_indices.push_back(_surfaces.size());
                guaranteed_index transform_index = _surface_transforms.size();
                _surface_transforms.push_back(std::move(transform));
                _surfaces.push_back(detector_surface(std::move(transform_index), std::move(typed_mask_range), std::move(volume_index), false));
            }
        }

        /** Const access method for the detector name */
        const std::string& name() const { return _name; }

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
    };

} // namespace detray