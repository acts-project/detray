/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Licenced under: Apache-2, see LICENSE file
 */
#pragma once

#include "core/intersection.hpp"

#include <array>
#include <cmath>
#include <climits>

namespace detray
{
    /** This is a simple 2-dimensional mask
     * 
     * It is defined by half length in local0 coord (h0) and local1 (h1), and can be checked
     * with a tolerance in t0 and t1.
     **/
    template <typename scalar_type>
    struct rectangle2
    {
        std::array<scalar_type,2> _h = 
            { std::numeric_limits<scalar_type>::infinity(),
              std::numeric_limits<scalar_type>::infinity() };

        /** Mask operation 
         * 
         * @tparam point2_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         * 
         * @param p the point to be checked
         * @param t0 is the tolerance in local 0
         * @param t1 is the tolerance in local 1
         * 
         * @return an intersection status e_inside / e_outside
         **/
        template <typename point2_type>
        intersection_status operator()(const point2_type& p, scalar_type t0=0., scalar_type t1=0.) const{
            return ( std::abs(p[0]) <= _h[0]+t0 and std::abs(p[1]) <= _h[1]+t1) ? e_inside : e_outside;
        }


        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const std::array<scalar_type, 2>& rhs){
            return (std::abs(_h[0]-rhs[0]) < std::numeric_limits<scalar_type>::epsilon()
                and std::abs(_h[1]-rhs[1]) < std::numeric_limits<scalar_type>::epsilon());
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const rectangle2<scalar_type>& rhs){
            return operator==(rhs._h);
        }

    };
    
} // namespace detray
