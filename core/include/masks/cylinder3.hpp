/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "utils/containers.hpp"
#include "tools/cylinder_intersector.hpp"

#include <cmath>
#include <climits>

namespace detray
{
    /** This is a simple mask for a full cylinder
     * 
     * It is defined by r and the half length.
     **/
    template <typename scalar_type, typename intersector_type = detray::cylinder_intersector>
    struct cylinder3
    {
        darray<scalar_type, 2> _v =
            {std::numeric_limits<scalar_type>::infinity(),
             std::numeric_limits<scalar_type>::infinity()};

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
                                       scalar_type t0 = std::numeric_limits<scalar_type>::epsilon(),
                                       scalar_type t1 = std::numeric_limits<scalar_type>::epsilon()) const
        {
            scalar_type r = getter::perp(p);
            if (std::abs(r - _v[0]) >= t0 + 5 * std::numeric_limits<scalar_type>::epsilon())
            {
                return e_missed;
            }
            return (std::abs(p[2]) <= _v[1] + t1) ? e_inside : e_outside;
        }

        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const darray<scalar_type, 2> &rhs)
        {
            return (std::abs(_v[0] - rhs[0]) < std::numeric_limits<scalar_type>::epsilon() and std::abs(_v[1] - rhs[1]) < std::numeric_limits<scalar_type>::epsilon());
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const cylinder3<scalar_type> &rhs)
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

        /** Return an associated intersector type */
        intersector_type intersector() { return intersector_type{}; };

    };

} // namespace detray
