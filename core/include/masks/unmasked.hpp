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

        using mask_tolerance = bool;

        static constexpr mask_tolerance within_epsilon = true;

        /** Mask operation 
         * 
         * @tparam point_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         *
         * the parameters are ignored
         * 
         * @return a bool that is ture if inside
         **/
        template <typename local_type>
        intersection_status is_inside(const typename local_type::point2 & /*ignored*/,
                                      const mask_tolerance& t = within_epsilon) const
        {
            return t ? e_hit : e_missed;
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
        template <typename local_type>
        bool is_inside(const typename local_type::point2 & /*ignored*/, scalar_type /*ignored*/) const
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
        template <typename local_type>
        bool is_inside(const typename local_type::point2 & /*ignored*/, scalar_type /*ignored*/, scalar_type /*ignored*/) const
        {
            return true;
        }

    };

} // namespace detray
