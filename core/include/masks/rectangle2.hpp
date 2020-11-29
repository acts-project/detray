#pragma once

#include "core/intersection.hpp"

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
        scalar_type _h0 = std::numeric_limits<scalar_type>::infinity();
        scalar_type _h1 = std::numeric_limits<scalar_type>::infinity();

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
            return ( std::abs(p[0]) <= _h0+t0 and std::abs(p[1]) <= _h1+t1) ? e_inside : e_outside;
        }
    };
    
} // namespace detray
