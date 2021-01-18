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

        /** Intersection method for cylindrical surfaces
         * 
         * @tparam transform_type The surface transform type to be intersected
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
         * 
         * @return the intersection with optional parameters
         **/
        template <typename transform_type, typename local_type, typename mask_type>
        intersection<scalar, typename transform_type::transform3::point3, typename local_type::point2>
        intersect(const transform_type &trf,
                  const track<transform_type> &track,
                  const local_type &local,
                  const mask_type &mask) const
        {
            return intersect(trf, track.pos, track.dir, local, mask, track.overstep_tolerance);
        }

        /** Intersection method for cylindrical surfaces
         * 
         * @tparam transform_type The transform type of the surface  to be intersected
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
         * 
         * @return the intersection with optional parameters
         **/
        template <typename transform_type, typename local_type, typename mask_type>
        intersection<scalar, typename transform_type::point3, typename local_type::point2>
        intersect(const transform_type &trf,
                  const typename transform_type::point3 &ro,
                  const typename transform_type::vector3 &rd,
                  const local_type &local,
                  const mask_type &mask,
                  scalar overstep_tolerance = 0.) const
        {
            using intersection = intersection<scalar, typename transform_type::point3, typename local_type::point2>;

            scalar r = mask[0];

            const auto &m = trf.matrix();
            auto sz = getter::vector<3>(m, 0, 2);
            auto sc = getter::vector<3>(m, 0, 3);

            typename transform_type::vector3 pc_cross_sz = vector::cross(ro - sc, sz);
            typename transform_type::vector3 rd_cross_sz = vector::cross(rd, sz);
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
                    is.point3 = ro + is.path * rd;
                    is.point2 = local(trf, is.point3);
                    auto local3 = trf.point_to_local(is.point3);
                    is.status = mask.template is_inside<transform_type>(local3);
                    scalar rdr = getter::perp(local3 + 10 * std::numeric_limits<scalar>::epsilon() * rd);
                    is.direction = rdr > r ? e_along : e_opposite;
                    return is;
                }
            }
            return intersection{};
        }
    };

} // namespace detray
