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
/** This is a simple mask for single parameter bound mask
 *
 * @tparam kCheckIndex is the index of the position on which the mask is applied
 * @tparam intersector_t is a struct used for intersecting this cylinder
 * @tparam local_type is the default local frame definition type
 * @tparam links_type is an object where the mask can link to
 * @tparam kMaskContext is a unique mask identifier in a certain context
 *
 * @note  While the mask_context can change depending on the typed container
 * structure the mask_identifier is a const expression that determines the
 * mask type once for all.
 *
 **/
template <unsigned int kCheckIndex,
          typename local_t = __plugin::cartesian2<detray::scalar>,
          typename links_t = dindex,
          template <typename, std::size_t> class array_t = darray>
class single3 final
    : public mask_base<plane_intersector, local_t, links_t, array_t, 2> {
    public:
    using base_type =
        mask_base<plane_intersector, local_t, links_t, array_t, 2>;
    using base_type::base_type;
    using mask_values = typename base_type::mask_values;
    using links_type = typename base_type::links_type;
    using local_type = typename base_type::local_type;
    using intersector_type = typename base_type::intersector_type;
    using point3 = __plugin::point3<scalar>;

    /* Default constructor */
    single3() : base_type({std::numeric_limits<scalar>::infinity()}, {}) {}

    /** Construction from boundary values **/
    DETRAY_HOST_DEVICE
    single3(scalar x, scalar y, links_type links) : base_type({x, y}, links) {}

    /** Assignment operator from an array, convenience function
     *
     * @param rhs is the right hand side object
     **/
    DETRAY_HOST_DEVICE
    single3<kCheckIndex, local_type, links_type, array_t> &operator=(
        const mask_values &rhs) {
        this->_values = rhs;
        return (*this);
    }

    /** Mask operation
     *
     * @tparam inside_local_type is the global type for checking (ignored)
     *
     * @param p the point to be checked
     * @param t is the tolerance of the single parameter
     *
     * @return an intersection status e_inside / e_outside
     **/
    template <typename inside_local_t>
    DETRAY_HOST_DEVICE intersection::status is_inside(
        const point3 &p,
        const scalar t = std::numeric_limits<scalar>::epsilon()) const {
        return (this->_values[0] - t <= p[kCheckIndex] and
                p[kCheckIndex] <= this->_values[1] + t)
                   ? intersection::status::e_inside
                   : intersection::status::e_outside;
    }

    /** Transform to a string for output debugging */
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "single3";
        for (const auto &v : this->_values) {
            ss << ", " << v;
        }
        return ss.str();
    }
};

}  // namespace detray
