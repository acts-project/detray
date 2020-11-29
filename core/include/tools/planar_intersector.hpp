#pragma once

#include <cmath>
#include <climits>

#include "core/intersection.hpp"
#include "core/surface.hpp"
#include "masks/unmasked.hpp"

namespace detray
{

    /** This is an intersector struct for planar surface
     */
    struct planar_intersector
    {
        
         /** A representation that is not bound to a local frame
         */
        struct unbound
        {
            using point2 = int;

             /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame,
              * including the contextual transform into the local 3D frame
              * 
              * @tparam the type of the surface from which also point3 and context type can be deduced
              * 
              */
            template <typename surface_type>
            const std::optional<point2> operator()(const surface_type & /*ignored*/, 
                                                   const typename surface_type::transform3::point3 & /*ignored*/,
                                                   const typename surface_type::transform3::context & /*ignored*/) const
            {
                return std::nullopt;
            }

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             */
            template <typename point3_type>
            const std::optional<point2> operator()(const point3_type & /*ignored*/) const
            {
                return std::nullopt;
            }
        };

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
         * 
         * Non-contextual part:
         * @param local the local local frame 
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
            auto sm = s.transform().matrix(ctx);
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
                is._status = mask(is._point2);

                return is;
            }
            return intersection{};
        }
    };

} // namespace detray
