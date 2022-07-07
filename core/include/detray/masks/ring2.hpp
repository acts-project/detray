/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/mask_base.hpp"

// System include(s)
#include <climits>
#include <cmath>
#include <sstream>
#include <string>

namespace detray {
/** This is a simple 2-dimensional mask for a closed ring
 *
 * @tparam intersector_t is a struct used for intersecting this cylinder
 * @tparam local_type is the default local frame type
 * @tparam links_type is an object where the mask can link to
 * @tparam kMaskContext is a unique mask identifier in a certain context
 *
 * It is defined by the two radii _values[0] and  _values[1],
 * and can be checked with a tolerance in t[0] and t[1].
 *
 * @note  While the mask_context can change depending on the typed container
 * structure the mask_identifier is a const expression that determines the
 * mask type once for all.
 *
 **/
template <typename local_t = __plugin::cartesian2<detray::scalar>,
          typename links_t = dindex,
          template <typename, std::size_t> class array_t = darray>
class ring2 final
    : public mask_base<plane_intersector, local_t, links_t, array_t, 2> {
    public:
    using base_type =
        mask_base<plane_intersector, local_t, links_t, array_t, 2>;
    using base_type::base_type;
    using mask_values = typename base_type::mask_values;
    using links_type = typename base_type::links_type;
    using local_type = typename base_type::local_type;
    using intersector_type = typename base_type::intersector_type;
    using point2 = __plugin::point2<scalar>;
    using point3 = __plugin::point3<scalar>;

    /* Default constructor */
    ring2() : base_type({0., std::numeric_limits<scalar>::infinity()}, {}) {}

    /** Construction from boundary values
     *
     * @param r_low lower radial bound
     * @param r_high upper radial bound
     */
    DETRAY_HOST_DEVICE
    ring2(scalar r_low, scalar r_high, links_type links)
        : base_type({r_low, r_high}, links) {}

    /** Assignment operator from an array, convenience function
     *
     * @param rhs is the right hand side object
     **/
    DETRAY_HOST_DEVICE
    ring2<local_type, links_type, array_t> &operator=(const mask_values &rhs) {
        this->_values = rhs;
        return (*this);
    }

    /** Mask operation
     *
     * @param p the point to be checked
     * @param t is the tolerance in r
     *
     * @return an intersection status e_inside / e_outside
     **/
    template <typename cartesian_point_t>
    DETRAY_HOST_DEVICE intersection::status is_inside(
        const cartesian_point_t &p,
        const scalar t = std::numeric_limits<scalar>::epsilon()) const {
        const scalar r = getter::perp(p);
        return (r + t >= this->_values[0] and r <= this->_values[1] + t)
                   ? intersection::status::e_inside
                   : intersection::status::e_outside;
    }

    /** Transform to a string for output debugging */
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "ring2";
        for (const auto &v : this->_values) {
            ss << ", " << v;
        }
        return ss.str();
    }
};

}  // namespace detray
