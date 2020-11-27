#pragma once

#include <climits>

namespace detray
{
    template <typename scalar_type>
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
        bool operator()(const point2_type & /*ignored*/) const
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