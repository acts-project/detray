/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/mask_base.hpp"

// System include(s)
#include <string>

namespace detray {

template <typename local_t = __plugin::cartesian2<detray::scalar>,
          typename links_t = dindex,
          template <typename, std::size_t> class array_t = darray>
class unmasked final
    : public mask_base<plane_intersector, local_t, links_t, array_t, 1> {
    public:
    using base_type =
        mask_base<plane_intersector, local_t, links_t, array_t, 1>;
    using base_type::base_type;
    using mask_values = typename base_type::mask_values;
    using links_type = typename base_type::links_type;
    using local_type = typename base_type::local_type;
    using intersector_type = typename base_type::intersector_type;
    using point2 = __plugin::point2<scalar>;

    /* Default constructor */
    unmasked() = default;

    /** Construction from voume link */
    DETRAY_HOST_DEVICE
    unmasked(links_type links) : base_type({0}, links) {}

    /** Mask operation
     *
     * @tparam point_type is the type of the point to be checked w.r.t. to
     * the mask bounds
     *
     * the parameters are ignored
     *
     * @return true
     **/
    template <typename inside_local_t>
    DETRAY_HOST_DEVICE inline intersection::status is_inside(
        const point2& /*ignored*/, const scalar /*ignored*/) const {
        return intersection::status::e_inside;
    }

    /** Transform to a string for output debugging */
    std::string to_string() const { return "unmasked"; }
};

}  // namespace detray
