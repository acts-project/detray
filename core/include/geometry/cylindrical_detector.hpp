/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "core/surface.hpp"
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
        // Indexing
        using optional_index = int;
        using guaranteed_index = unsigned long;

        /** Nested volume struct */
        struct volume
        {
            dvector<guaranteed_index> portal_surface_indices = {};
            guaranteed_index volume_index = 0;
        };

        // Portal handling
        using portal_destinations = darray<optional_index, 2>;
        using portal_cylinder_mask = cylinder3<scalar, false, concentric_cylinder_intersector, portal_destinations>;
        using portal_cylinders = dvector<portal_cylinder_mask>;
        using portal_disc_mask = ring2<scalar, planar_intersector, portal_destinations>;
        using portal_discs = dvector<portal_disc_mask>;
        using portal_type_map = dmap<guaranteed_index, guaranteed_index>;
        using portal_range = darray<guaranteed_index, 2>;
        using portal_links = dtuple<guaranteed_index, portal_range>;
        using portal_surface = surface<transform_type, portal_links, guaranteed_index>;

        // Member & method section: for portals
    private:
        dtuple<portal_cylinders, portal_discs> _portals;
        dvector<portal_surface> _portal_surfaces;
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
            portal_range range = {index_start, index_start + portals.size()};
            guaranteed_index type = _portal_types.find(portal_type::mask_identifier)->first;
            // Record the portal index
            volume.portal_surface_indices.push_back(_portal_surfaces.size());
            guaranteed_index index = volume.volume_index;
            portal_links links = {type, range};
            _portal_surfaces.push_back(portal_surface(std::move(transform), std::move(links), std::move(index), false));
        }

        // Member & method section for volumes
    private:
        dvector<volume> _volumes;

    public:
        /** Add a volume to the detector 
         * 
         * @tparam volume the volume to be added to the detector
         **/
        void add_volume(const volume &volume)
        {
            _volumes.push_back(volume);
        }

    };

} // namespace detray