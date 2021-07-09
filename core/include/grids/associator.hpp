/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <vector>
#include <climits>

namespace detray
{

    using point2 = __plugin::point2;

    /** Struct that assigns the center of gravity to a rectangular bin */
    struct center_of_gravity_rectangle
    {
        /** Call operator to the struct, allows to chain several chain operators together
         * 
         * @param bin_contour The contour description of the bin -> target
         * @param surface_contour The contour description of the surface -> test
         * 
         * @note the bin_contour is asummed to be a rectangle
         * 
         * @return whether this should be associated
         */
        bool operator()(const std::vector<point2> &bin_contour, const std::vector<point2> &surface_contour)
        {
            // Check if centre of gravity is inside bin
            point2 cgs = {0., 0.};
            for (const auto &svtx : surface_contour)
            {
                cgs = cgs + svtx;
            }
            cgs = 1. / surface_contour.size() * cgs;
            scalar min_l0 = std::numeric_limits<scalar>::max();
            scalar max_l0 = -std::numeric_limits<scalar>::max();
            scalar min_l1 = std::numeric_limits<scalar>::max();
            scalar max_l1 = -std::numeric_limits<scalar>::max();
            for (const auto &b : bin_contour)
            {
                min_l0 = std::min(b[0], min_l0);
                max_l0 = std::max(b[0], max_l0);
                min_l1 = std::min(b[1], min_l1);
                max_l1 = std::max(b[1], max_l1);
            }

            if (cgs[0] >= min_l0 and cgs[0] < max_l0 and cgs[1] >= min_l1 and cgs[1] < max_l1)
            {
                return true;
            }

            return false;
        }
    };

    /** Check if center of mass is inside a generic polygon bin */
    struct center_of_gravity_generic
    {
        /** Call operator to the struct, allows to chain several chain operators together
        * 
        * @param bin_contour The contour description of the bin -> target
        * @param surface_contour The contour description of the surface -> test
        * 
        * @return whether this should be associated
         */
        bool operator()(const std::vector<point2> &bin_contour, const std::vector<point2> &surface_contour)
        {
            // Check if centre of gravity is inside bin
            point2 cgs = {0., 0.};
            for (const auto &svtx : surface_contour)
            {
                cgs = cgs + svtx;
            }
            cgs = 1. / surface_contour.size() * cgs;

            size_t i, j = 0;
            size_t num_points = bin_contour.size();

            bool inside = false;
            for (i = 0, j = num_points - 1; i < num_points; j = i++)
            {
                const auto &pi = bin_contour[i];
                const auto &pj = bin_contour[j];
                if ((((pi[1] <= cgs[1]) and (cgs[1] < pj[1])) or
                     ((pj[1] <= cgs[1]) and (cgs[1] < pi[1]))) and
                    (cgs[0] < (pj[0] - pi[0]) * (cgs[1] - pi[1]) / (pj[1] - pi[1]) + pi[0]))
                    inside = !inside;
            }
            return inside;
        }
    };

     /** Check if the egdes of the bin and surface contour overlap */
    struct edges_intersect_generic
    {

        /** Call operator to the struct, allows to chain several chain operators together
       * 
       * @param bin_contour The contour description of the bin -> target
       * @param surface_contour The contour description of the surface -> test
       * 
       * @return whether this should be associated
       */
        bool operator()(const std::vector<point2> &bin_contour, const std::vector<point2> &surface_contour)
        {

            auto intersect = [](const point2 &pi, const point2 &pj, const point2 &pk, const point2 &pl) -> bool
            {
                scalar d = (pj[0] - pi[0]) * (pl[1] - pk[1]) - (pj[1] - pi[1]) * (pl[0] - pk[0]);

                if (d != 0.)
                {
                    double r = ((pi[1] - pk[1]) * (pl[0] - pk[0]) - (pi[0] - pk[0]) * (pl[1] - pk[1])) / d;
                    double s = ((pi[1] - pk[1]) * (pj[0] - pi[0]) - (pi[0] - pk[0]) * (pj[1] - pi[1])) / d;
                    if (r >= 0. and r <= 1. and s >= 0. and s <= 1.)
                    {
                        return true;
                    }
                }
                return false;
            };

            // Loop over bin_contour
            for (size_t j = 1; j <= bin_contour.size(); ++j)
            {
                size_t i = j - 1;
                size_t jc = (j == bin_contour.size()) ? 0 : j;
                const auto &pi = bin_contour[i];
                const auto &pj = bin_contour[jc];
                // Loop over surface_contour
                for (size_t k = 1; k <= surface_contour.size(); ++k)
                {
                    size_t l = k - 1;
                    size_t kc = (k == surface_contour.size()) ? 0 : k;
                    const auto &pl = surface_contour[l];
                    const auto &pk = surface_contour[kc];
                    if (intersect(pi, pj, pk, pl))
                    {
                        return true;
                    }
                }
            }
            return false;
        }
    };

} // namespace detray