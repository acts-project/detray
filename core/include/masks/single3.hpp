/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "utils/containers.hpp"

#include <cmath>
#include <climits>
#include <iostream>

namespace detray
{
    /** This is a simple mask for single parameter bound mask
     * 
     * It is defined by r and the half length.
     **/
    template <typename scalar_type, unsigned int kPAR>
    struct single3
    {
        darray<scalar_type, 1> _v =
            {std::numeric_limits<scalar_type>::infinity()};

        /** Mask operation 
         * 
         * @tparam point3_type is the type of the point to be checked w.r.t. to
         * the mask bounds, it's assumed to be within the cylinder 3D frame
         * 
         * @param p the point to be checked
         * @param t0 is the tolerance in local 0 (radius)
         * @param t1 is the tolerance in local 1 (z)
         * 
         * @return an intersection status e_inside / e_outside
         **/
        template <typename point3_type>
        intersection_status operator()(const point3_type &p,
                                       scalar_type t = std::numeric_limits<scalar_type>::epsilon()) const
        {     
            return (std::abs(p[kPAR]) <= _v[0] + t) ? e_inside : e_outside;
        }

        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const darray<scalar_type, 2> &rhs)
        {
            return (std::abs(_v[0] - rhs[0]) < std::numeric_limits<scalar_type>::epsilon());
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const single3<scalar_type, kPAR> &rhs)
        {
            return operator==(rhs._v);
        }

        /** Access operator - non-const
         * @return the reference to the member variable
         */
        scalar_type &operator[](unsigned int value_index)
        {
            return _v[value_index];
        }

        /** Access operator - non-const
         * @return a copy of the member variable
         */
        scalar_type operator[](unsigned int value_index) const
        {
            return _v[value_index];
        }
    };

} // namespace detray
