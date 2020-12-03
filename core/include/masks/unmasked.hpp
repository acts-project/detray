/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"

#include <optional>

namespace detray
{
    template <typename scalar_type = std::nullopt_t>
    struct unmasked
    {
        
        /** Mask operation 
         * 
         * @tparam point_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         *
         * the parameters are ignored
         * 
         * @return a bool that is ture if inside
         **/
        template <typename point_type>
        intersection_status operator()(const point_type & /*ignored*/) const
        {
            return e_hit;
        }

        /** Mask operation 
         * 
         * @tparam point_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         *
         * the parameters are ignored
         * 
         * @return a bool that is ture if inside
         **/
        template <typename point_type>
        bool operator()(const point_type & /*ignored*/, scalar_type /*ignored*/) const
        {
            return true;
        }

        /** Mask operation 
         * 
         * @tparam point_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         *
         * the parameters are ignored
         * 
         * @return a bool that is ture if inside
         **/
        template <typename point_type>
        bool operator()(const point_type & /*ignored*/, scalar_type /*ignored*/, scalar_type /*ignored*/) const
        {
            return true;
        }
    };

} // namespace detray