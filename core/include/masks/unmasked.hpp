#pragma once

#include "core/intersection.hpp"

#include <climits>
#include <optional>

namespace detray
{
    template <typename scalar_type = std::nullopt_t>
    struct unmasked
    {
        
        /** Mask operation 
         * 
         * @tparam point2_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         *
         * the parameters are ignored
         * 
         * @return a bool that is ture if inside
         **/
        template <typename point2_type>
        intersection_status operator()(const point2_type & /*ignored*/) const
        {
            return e_hit;
        }

        /** Mask operation 
         * 
         * @tparam point2_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         *
         * the parameters are ignored
         * 
         * @return a bool that is ture if inside
         **/
        template <typename point2_type>
        bool operator()(const point2_type & /*ignored*/, scalar_type /*ignored*/) const
        {
            return true;
        }

        /** Mask operation 
         * 
         * @tparam point2_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         *
         * the parameters are ignored
         * 
         * @return a bool that is ture if inside
         **/
        template <typename point2_type>
        bool operator()(const point2_type & /*ignored*/, scalar_type /*ignored*/, scalar_type /*ignored*/) const
        {
            return true;
        }
    };

} // namespace detray