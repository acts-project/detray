/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinates.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/mask_base.hpp"

// System include(s).
#include <climits>
#include <cmath>
#include <sstream>
#include <string>

namespace detray {

/** This is a simple 2-dimensional mask for a regular rectangle
 *
 * @tparam intersector_t is a struct used for intersecting this cylinder
 * @tparam local_type is the default local frame definition type
 * @tparam links_type is an object where the mask can link to
 * @tparam kMaskContext is a unique mask identifier in a certain context
 *
 * It is defined by half length in local0 coordinates _values[0] and _values[1],
 * and can be checked with a tolerance in t[0] and t[1].
 *
 * @note  While the mask_context can change depending on the typed container
 * structure the mask_identifier is a const expression that determines the
 * mask type once for all.
 *
 **/

template <typename transform3_t = __plugin::transform3<scalar>,
          template <class> typename intersector_t = plane_intersector,
          template <class> typename local_t = cartesian2,
          typename links_t = dindex,
          template <typename, std::size_t> class array_t = darray>
class rectangle2 final : public mask_base<transform3_t, intersector_t, local_t,
                                          links_t, array_t, 2> {
    public:
    using base_type =
        mask_base<transform3_t, intersector_t, local_t, links_t, array_t, 2>;
    using base_type::base_type;
    using mask_values = typename base_type::mask_values;
    using links_type = typename base_type::links_type;
    using local_type = typename base_type::local_type;
    using intersector_type = typename base_type::intersector_type;
    using point2 = typename transform3_t::point2;

    /* Default constructor */
    rectangle2()
        : base_type({std::numeric_limits<scalar>::infinity(),
                     std::numeric_limits<scalar>::infinity()},
                    {}) {}

    /** Construction from boundary values
     *
     * @param half_length_0 half length in loc0
     * @param half_length_1 half length in loc1
     */
    DETRAY_HOST_DEVICE
    rectangle2(scalar half_length_0, scalar half_length_1, links_type links)
        : base_type({half_length_0, half_length_1}, links) {}

    /** Mask operation
     *
     * @tparam inside_local_t is the local type for inside checking
     *
     * @param p the point to be checked
     * @param t is the tolerance tuple in (l0, l1)
     *
     * @return an intersection status e_inside / e_outside
     **/
    template <typename inside_local_t>
    DETRAY_HOST_DEVICE intersection::status is_inside(
        const point2 &p,
        const scalar t = std::numeric_limits<scalar>::epsilon()) const {
        return (std::abs(p[0]) <= this->_values[0] + t and
                std::abs(p[1]) <= this->_values[1] + t)
                   ? intersection::status::e_inside
                   : intersection::status::e_outside;
    }

    /** Transform to a string for output debugging */
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "rectangle2";
        for (const auto &v : this->_values) {
            ss << ", " << v;
        }
        return ss.str();
    }
};

}  // namespace detray
