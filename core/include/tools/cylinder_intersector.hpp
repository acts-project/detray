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
        
        using transform3 = __plugin::transform3;
        using point3 = __plugin::point3;
        using vector3 = __plugin::vector3;
        using point2 = __plugin::point2;
        using cylindrical2 = __plugin::cylindrical2;

        /** Intersection method for cylindrical surfaces
         * 
         * @tparam track_type The type of the track carrying the context object
         * @tparam local_type The local frame type to be intersected
         * @tparam mask_type The mask type applied to the local frame
         * 
         * Contextual part:
         * @param trf the transform of the surface to be intersected
         * @param track the track information
         * @param local to the local local frame 
         * 
         * Non-contextual part:
         * @param mask the local mask 
         * @param tolerance is the mask specific tolerance
         * 
         * @return the intersection with optional parameters
         **/
        template <typename track_type, typename local_type, typename mask_type>
        intersection
        intersect(const transform3 &trf,
                  const track_type &track,
                  const local_type &local,
                  const mask_type &mask,
                  const typename mask_type::mask_tolerance &tolerance = mask_type::within_epsilon) const
        {
            return intersect(trf, track.pos, track.dir, local, mask, tolerance, track.overstep_tolerance);
        }

        /** Intersection method for cylindrical surfaces
         * 
         * @tparam local_type The local frame type to be intersected
         * @tparam mask_type The mask type applied to the local frame
         * 
         * Contextual part:
         * @param trf the transform of the surface to be intersected
         * @param ro the origin of the ray
         * @param rd the direction of the ray
         * @param local to the local local frame 
         * 
         * Non-contextual part:
         * @param mask the local mask 
         * @param tolerance is the mask specific tolerance
         * @param overstep_tolerance  is the stepping specific tolerance
         * 
         * @return the intersection with optional parameters
         **/
        template <typename local_type, typename mask_type>
        intersection
        intersect(const transform3 &trf,
                  const point3 &ro,
                  const vector3 &rd,
                  const local_type &local,
                  const mask_type &mask,
                  const typename mask_type::mask_tolerance &tolerance = mask_type::within_epsilon,
                  scalar overstep_tolerance = 0.) const
        {

            scalar r = mask[0];

            const auto &m = trf.matrix();
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
                    is.path = t;
                    is.p3 = ro + is.path * rd;
                    is.p2 = local(trf, is.p3);
                    auto local3 = trf.point_to_local(is.p3);
                    is.status = mask.template is_inside<cylindrical2>(local3, tolerance);
                    scalar rdr = getter::perp(local3 + 10 * std::numeric_limits<scalar>::epsilon() * rd);
                    is.direction = rdr > r ? e_along : e_opposite;
                    return is;
                }
            }
            return intersection{};
        }
    };

} // namespace detray
