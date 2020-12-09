/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <cmath>
#include <climits>
#include <optional>

#include "core/intersection.hpp"
#include "core/surface.hpp"
#include "masks/unmasked.hpp"
#include "utils/unbound.hpp"

namespace detray
{

    /** This is an intersector struct for planar surface
     */
    struct planar_intersector
    {
        
        /** Intersection method for planar surfaces
         * 
         * @tparam surface_type The surface type to be intersected
         * @tparam local_type The local frame type to be intersected
         * @tparam mask_type The mask type applied to the local frame
         * 
         * Contextual part:
         * @param s the surface to be intersected
         * @param ro the origin of the ray
         * @param rd the direction of the ray
         * @param ctx the context of the call
         * @param local to the local local frame 
         * 
         * Non-contextual part:
         * @param mask the local mask 
         * 
         * @return the intersection with optional parameters
         **/
        template <typename surface_type, typename local_type = unbound, typename mask_type = unmasked<>>
        intersection<scalar, typename surface_type::transform3::point3, typename local_type::point2>
        intersect(const surface_type &s,
                  const typename surface_type::transform3::point3 &ro,
                  const typename surface_type::transform3::vector3 &rd,
                  const typename surface_type::transform3::context &ctx,
                  const local_type &local = local_type(),
                  const mask_type &mask = mask_type()) const
        {

            using point3 = typename surface_type::transform3::point3;
            using point2 = typename local_type::point2;
            using intersection = intersection<scalar, point3, point2>;

            // Retrieve the surface normal & translation (context resolved)
            const auto& sm = s.transform().matrix(ctx);
            auto sn = getter::block<3, 1>(sm, 0, 2);
            auto st = getter::block<3, 1>(sm, 0, 3);

            // Intersection code
            scalar denom = vector::dot(rd, sn);
            if (denom != 0.0)
            {
                intersection is;
                is._path = vector::dot(sn, (st - ro)) / (denom);
                is._point3 = ro + is._path * rd;
                is._point2 = local(s, is._point3, ctx);
                is._status = mask(is._point2.value_or(point2()));
                is._direction = denom > 0 ? e_along : e_opposite;

                return is;
            }
            return intersection{};
        }
    };

} // namespace detray
