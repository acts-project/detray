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

    template <typename transform_type, typename source_link = bool>
    struct cylindrical_detector
    {
        using optional_index = int;
        using guaranteed_index = unsigned long;

        using portal_destinations = darray<int, 2>;
        using portal_cylinder_mask = cylinder3<scalar, false, concentric_cylinder_intersector, portal_destinations>;
        using portal_cylinders = dvector<portal_cylinder_mask>;
        using portal_disc_mask = ring2<scalar, planar_intersector, portal_destinations>;
        using portal_discs = dvector<portal_disc_mask>;
        using portal_type_map = dmap<unsigned int, guaranteed_index>;
        using portal_range = std::array<guaranteed_index, 2>;
        using portal_links = dtuple<guaranteed_index, portal_range>;
        using portal_surface = surface<transform_type, portal_links>;
        using portal_surface_range = std::array<guaranteed_index, 2>;

        dtuple<portal_cylinders, portal_discs> _portals;
        dvector<portal_surface> _portal_surfaces;
        portal_type_map _portal_types = {{portal_cylinder_mask::mask_identifier, 0}, {portal_disc_mask::mask_identifier, 1}};

        struct volume
        {
            portal_surface_range portal_surface_range = {0, 0};
            guaranteed_index volume_index = 0;
        };

        //using mask_index = std::array<int, 2>;

        //dvector<volume> _volumes;
        //dvector<surface> _surfaces;

        /** Method to add a list of portals to the portal tuple
         * 
         * @tparam portal_type the type of the portal mask
         * 
         * @param portals the vector of portal masks
         * 
         * It @returns an indexed portal group
         **/
        template <typename portal_type>
        portal_links add_portals(const dvector<portal_type> &portals)
        {
            auto &group = std::get<dvector<portal_type>>(_portals);
            guaranteed_index index_start = group.size();
            portal_range range = {index_start, index_start + portals.size()};
            group.insert(group.end(), portals.begin(), portals.end());
            guaranteed_index type = _portal_types.find(portal_type::mask_identifier)->first;
            return {type, range};
        }

        guaranteed_index add_portal_surface(portal_surface&& psurface) 
        {
          guaranteed_index this_index = _portal_surfaces.size();
          _portal_surfaces.push_back(std::move(psurface));
          return this_index;
        } 

    };

} // namespace detray