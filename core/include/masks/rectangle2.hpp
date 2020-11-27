#pragma once

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
        scalar_type h0 = std::numeric_limits<scalar_type>::infinity();
        scalar_type h1 = std::numeric_limits<scalar_type>::infinity();

        /** Mask operation 
         * 
         * @tparam point2_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         * 
         * @param p the point to be checked
         * @param t0 is the tolerance in local 0
         * @param t1 is the tolerance in local 1
         * 
         * @return a bool that is ture if inside
         **/
        template <typename point2_type>
        bool operator()(const point2_type& p, scalar_type t0=0., scalar_type t1=0.) const{
            return std::abs(p[0]) <= h0+t0 and std::abs(p[1]) <= h1+t1;
        }
    };
    
} // namespace detray
