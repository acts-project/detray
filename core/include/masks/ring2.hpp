/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"

#include <array>
#include <cmath>
#include <climits>

namespace detray
{
    /** This is a simple 2-dimensional mask for a closed ring
     * 
     * It is defined by the two radii _r[0] and  _r[1], 
     * and can be checked with a tolerance in t0 and t1.
     **/
    template <typename scalar_type>
    struct ring2
    {
        std::array<scalar_type, 2> _r =
            {0.,
             std::numeric_limits<scalar_type>::infinity()};

        /** Mask operation 
         * 
         * @tparam point2_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         * 
         * @param p the point to be checked
         * @param t0 is the tolerance in local 0
         * @param t1 is the tolerance in local 1 and is ignored
         * 
         * @return an intersection status e_inside / e_outside
         **/
        template <typename point2_type>
        intersection_status operator()(const point2_type &p,
                                       scalar_type t0 = std::numeric_limits<scalar_type>::epsilon(),
                                       scalar_type t1 = std::numeric_limits<scalar_type>::epsilon()) const
        {
            scalar_type r = vector::perp(p - point2_type{t0, t1});
            return (r + t0 >= _r[0] and r <= _r[1] + t0) ? e_inside : e_outside;
        }

        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const std::array<scalar_type, 2> &rhs)
        {
            return (std::abs(_r[0] - rhs[0]) < std::numeric_limits<scalar_type>::epsilon() and std::abs(_r[1] - rhs[1]) < std::numeric_limits<scalar_type>::epsilon());
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const ring2<scalar_type> &rhs)
        {
            return operator==(rhs._r);
        }

        /** Access operator - non-const
         * @return the reference to the member variable
         */
        scalar_type& operator[](unsigned int value_index) {
            return _r[value_index];
        }

        /** Access operator - non-const
         * @return a copy of the member variable
         */
        scalar_type operator[](unsigned int value_index) const {
            return _r[value_index];
        }

    };

} // namespace detray
