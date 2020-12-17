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
#include "core/track.hpp"
#include "utils/quadratic_equation.hpp"
#include "utils/unbound.hpp"
#include "masks/unmasked.hpp"

namespace detray
{

    /** This is an intersector struct for generic cylinder surface
     */
    struct cylinder_intersector
    {

        /** Intersection method for cylindrical surfaces
         * 
         * @tparam surface_type The surface type to be intersected
         * @tparam local_type The local frame type to be intersected
         * @tparam mask_type The mask type applied to the local frame
         * 
         * Contextual part:
         * @param s the surface to be intersected
         * @param t the track information
         * @param local to the local local frame 
         * 
         * Non-contextual part:
         * @param mask the local mask 
         * 
         * @return the intersection with optional parameters
         **/
        template <typename surface_type, typename local_type, typename mask_type>
        intersection<scalar, typename surface_type::transform3::point3, typename local_type::point2>
        intersect(const surface_type &s,
                  const track<typename surface_type::transform3> &t,
                  const local_type &local,
                  const mask_type &mask) const
        {
            return intersect(s, t.pos, t.dir, t.ctx, local, mask);
        }

        /** Intersection method for cylindrical surfaces
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
        template <typename surface_type, typename local_type, typename mask_type>
        intersection<scalar, typename surface_type::transform3::point3, typename local_type::point2>
        intersect(const surface_type &s,
                  const typename surface_type::transform3::point3 &ro,
                  const typename surface_type::transform3::vector3 &rd,
                  const typename surface_type::transform3::context &ctx,
                  const local_type &local,
                  const mask_type &mask) const
        {
            using point3 = typename surface_type::transform3::point3;
            using vector3 = typename surface_type::transform3::vector3;
            using point2 = typename local_type::point2;
            using intersection = intersection<scalar, point3, point2>;

            scalar r = mask[0];

            const auto &m = s.transform().matrix(ctx);
            auto sz = getter::vector<3>(m, 0, 2);
            auto sc = getter::vector<3>(m, 0, 3);

            vector3 pc_cross_sz = vector::cross(ro - sc, sz);
            vector3 rd_cross_sz = vector::cross(rd, sz);
            scalar a = vector::dot(rd_cross_sz, rd_cross_sz);
            double b = 2. * vector::dot(rd_cross_sz, pc_cross_sz);
            double c = vector::dot(pc_cross_sz, pc_cross_sz) - (r * r);

            quadratic_equation<scalar> qe = {a, b, c};
            auto qe_solution = qe();

            if (std::get<0>(qe_solution) > 0)
            {
                auto t01 = std::get<1>(qe_solution);
                scalar t = (t01[0] > 0.) ? t01[0] : t01[1];
                if (t > 0)
                {
                    intersection is;
                    is._path = t;
                    is._point3 = ro + is._path * rd;
                    is._point2 = local(s, is._point3, ctx);
                    auto local3 = s.transform().point_to_local(is._point3, ctx);
                    is._status = mask(local3);
                    scalar rdr = getter::perp(local3 + 10 * std::numeric_limits<scalar>::epsilon() * rd);
                    is._direction = rdr > r ? e_along : e_opposite;
                    return is;
                }
            }
            return intersection{};
        }
    };

} // namespace detray
