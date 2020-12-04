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
         * @tparam surface_type The surface type to be intersected
         * @tparam local_type The local frame type to be intersected
         * @tparam mask_type The mask type applied to the local frame
         * 
         * Contextual part:
         * @param s the surface to be intersected
         * @param r the radius of the surface
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
        intersect(const surface_type &s, scalar r,
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

            // Two points on the line
            const auto &l0 = ro;
            const auto l1 = ro + rd;

            // swap coorinates x/y for numerical stability
            bool swap_x_y = std::abs(rd[0]) < 1e-3;

            unsigned int _x = swap_x_y ? 1 : 0;
            unsigned int _y = swap_x_y ? 0 : 1;

            scalar k = (l0[_y] - l1[_y]) / (l0[_x] - l1[_x]);
            scalar d = l1[_y] - k * l1[_x];

            quadratic_equation<scalar> qe = {(1 + k * k), 2 * k, d * d - r * r};
            auto qe_solution = qe();

            if (std::get<0>(qe_solution) > 0)
            {

                darray<point3, 2> candidates;
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

                // chosse the index, take the smaller positive one
                int cindex = (t01[0] < t01[1] and t01[0] > 0.) ? 0 : (t01[0] < 0. and t01[1] > 0. ? 1 : -1);

                if (cindex > 0)
                {
                    intersection is;
                    is._point3 = candidates[cindex];
                    is._path = t01[cindex];
                    is._point2 = std::nullopt;
                    is._status = mask(is._point3);
                    return is;
                }
            }
            return intersection{};
        }
    };

} // namespace detray
