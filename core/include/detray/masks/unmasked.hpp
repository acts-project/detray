/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <optional>
#include <string>

#include "detray/core/intersection.hpp"
#include "detray/utils/unbound.hpp"

namespace detray {
template<typename mask_local_type = unbound>
struct unmasked {

    using mask_tolerance = bool;
    using local_type = mask_local_type;

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
    inline intersection_status is_inside(
        const point2 & /*ignored*/,
        const mask_tolerance &t = within_epsilon) const {
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
    inline bool is_inside(const point2 & /*ignored*/,
                          scalar /*ignored*/) const {
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
    inline bool is_inside(const point2 & /*ignored*/, scalar /*ignored*/,
                          scalar /*ignored*/) const {
        return true;
    }

    /** Transform to a string for output debugging */
    std::string to_string() const { return "unmasked"; }
};

}  // namespace detray
