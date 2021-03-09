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
    /** This is an intersector struct for a concetric cylinder surface
     */
    struct concentric_cylinder_intersector
    {

        /** Intersection method for cylindrical surfaces
         * 
         * @tparam transform_type The transform type of the surface to be intersected
         * @tparam local_type The local frame type to be intersected
         * @tparam mask_type The mask type applied to the local frame
         * 
         * Contextual part:
         * @param trf the transform of the surface surface to be intersected @note is ignored 
         * @param track the track information
         * @param local to the local local frame 
         * 
         * Non-contextual part:
         * @param mask the local mask 
         * @param tolerance is the mask specific tolerance
         * 
         * @return the intersection with optional parameters
         **/
        template <typename transform_type, typename local_type, typename mask_type>
        intersection<scalar, typename transform_type::point3, typename local_type::point2>
        intersect(const transform_type &trf,
                  const track<transform_type> &track,
                  const local_type &local,
                  const mask_type &mask,
                  const typename mask_type::mask_tolerance &tolerance = mask_type::within_epsilon) const
        {
            return intersect(trf, track.pos, track.dir, local, mask, track.overstep_tolerance);
        }

        /** Intersection method for cylindrical surfaces
         * 
         * @tparam transform_type The surface transform type to be intersected
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
        template <typename transform_type, typename local_type, typename mask_type>
        intersection<scalar, typename transform_type::point3, typename local_type::point2>
        intersect(const transform_type & /*trf*/,
                  const typename transform_type::point3 &ro,
                  const typename transform_type::vector3 &rd,
                  const local_type &local,
                  const mask_type &mask,
                  const typename mask_type::mask_tolerance &tolerance = mask_type::within_epsilon,
                  scalar overstep_tolerance = 0.) const
        {
            using intersection = intersection<scalar, typename transform_type::point3, typename local_type::point2>;

            scalar r = mask[0];

            // Two points on the line, thes are in the cylinder frame
            const auto &l0 = ro;
            const auto l1 = typename transform_type::point3(ro + rd);

            // swap coorinates x/y for numerical stability
            bool swap_x_y = std::abs(rd[0]) < 1e-3;

            unsigned int _x = swap_x_y ? 1 : 0;
            unsigned int _y = swap_x_y ? 0 : 1;
            scalar k = (l0[_y] - l1[_y]) / (l0[_x] - l1[_x]);
            scalar d = l1[_y] - k * l1[_x];

            quadratic_equation<scalar> qe = {(1 + k * k), 2 * k * d, d * d - r * r};
            auto qe_solution = qe();

            if (std::get<0>(qe_solution) > 0)
            {
                darray<typename transform_type::point3, 2> candidates;
                auto u01 = std::get<1>(qe_solution);
                darray<scalar, 2> t01 = {0., 0.};

                candidates[0][_x] = u01[0];
                candidates[0][_y] = k * u01[0] + d;
                t01[0] = (candidates[0][_x] - ro[_x]) / rd[_x];
                candidates[0][2] = ro[2] + t01[0] * rd[2];

                candidates[1][_x] = u01[1];
                candidates[1][_y] = k * u01[1] + d;
                t01[1] = (candidates[1][_x] - ro[_x]) / rd[_x];
                candidates[1][2] = ro[2] + t01[1] * rd[2];

                // Chose the index, take the smaller positive one
                int cindex = (t01[0] < t01[1] and t01[0] > overstep_tolerance) ? 0
                                                                               : (t01[0] < overstep_tolerance and t01[1] > 0. ? 1 : 0);
                if (t01[0] > overstep_tolerance or t01[1] > overstep_tolerance)
                {
                    intersection is;
                    is.point3 = candidates[cindex];
                    is.path = t01[cindex];

                    is.point2 = typename local_type::point2{r * getter::phi(is.point3), is.point3[2]};
                    is.status = mask.template is_inside<transform_type>(is.point3);
                    scalar rdir = getter::perp(is.point3 + 0.1 * rd);
                    is.direction = rdir > r ? e_along : e_opposite;
                    return is;
                }
            }
            return intersection{};
        }
    };

} // namespace detray
