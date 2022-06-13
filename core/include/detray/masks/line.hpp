/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <climits>
#include <cmath>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/ray_line_intersector.hpp"
#include "detray/masks/mask_base.hpp"

namespace detray {

/** This is a simple mask for a line which is defined with a line length and its
 * radial scope
 *
 **/
template <typename intersector_t = ray_line_intersector,
          typename local_t = __plugin::cartesian2<detray::scalar>,
          typename links_t = dindex,
          template <typename, std::size_t> class array_t = darray>
class line final
    : public mask_base<intersector_t, local_t, links_t, array_t, 2> {
    public:
    using base_type = mask_base<intersector_t, local_t, links_t, array_t, 2>;
    using base_type::base_type;
    using mask_tolerance = array_t<scalar, 2>;
    using mask_values = typename base_type::mask_values;
    using links_type = typename base_type::links_type;
    using local_type = typename base_type::local_type;
    using intersector_type = typename base_type::intersector_type;
    using point2 = __plugin::point2<scalar>;

    static constexpr mask_tolerance within_epsilon = {
        std::numeric_limits<scalar>::epsilon(),
        std::numeric_limits<scalar>::epsilon()};

    /* Default constructor */
    line() : base_type({0, std::numeric_limits<scalar>::infinity()}, {}) {}

    /** Construction from boundary values
     *
     * @param half_length is the half length of line
     * @param scope is the radial scope length of line detector
     */
    DETRAY_HOST_DEVICE
    line(scalar scope, scalar half_length, links_type links)
        : base_type({scope, half_length}, links) {}

    /** Assignment operator from an array, convenience function
     *
     * @param rhs is the right hand side object
     **/
    DETRAY_HOST_DEVICE
    line<intersector_t, local_type, links_type, array_t> &operator=(
        const array_t<scalar, 2> &rhs) {
        this->_values = rhs;
        return (*this);
    }

    /** Mask operation
     *
     * @tparam inside_local_t is the local type for inside checking
     *
     * @param p the point to be checked
     * @param t is the tolerance in r
     *
     * @return an intersection status e_inside / e_outside
     **/
    template <typename inside_local_t>
    DETRAY_HOST_DEVICE intersection::status is_inside(
        const point2 &p, const mask_tolerance t = within_epsilon) const {

        // if the point is whithin the scope, return e_inside
        return (std::abs(p[0]) <= this->_values[0] + t[0] &&
                std::abs(p[1]) <= this->_values[1] + t[1])
                   ? intersection::status::e_inside
                   : intersection::status::e_outside;
    }
};

}  // namespace detray