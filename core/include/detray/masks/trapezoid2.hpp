/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <climits>
#include <cmath>
#include <sstream>
#include <string>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/ray_plane_intersector.hpp"
#include "detray/masks/mask_base.hpp"

namespace detray {
/** This is a simple 2-dimensional mask for a regular trapezoid
 *
 * @tparam intersector_t is a struct used for intersecting this cylinder
 * @tparam local_type is the default local frame definition type
 * @tparam links_type is an object where the mask can link to
 * @tparam kMaskContext is a unique mask identifier in a certain context
 *
 * It is defined by half lengths in local0 coordinate _values[0] and _values[1]
 * at -/+ half length in the local1 coordinate _values[2]. _values[3] contains
 * the precomputed value of 1 / (2 * _values[2]), which avoids
 * excessive floating point divisions.
 *
 * @note  While the mask_context can change depending on the typed container
 * structure the mask_identifier is a const expression that determines the
 * mask type once for all.
 *
 **/
template <typename intersector_t = ray_plane_intersector,
          typename local_t = __plugin::cartesian2<detray::scalar>,
          typename links_t = dindex,
          template <typename, std::size_t> class array_t = darray>
class trapezoid2 final
    : public mask_base<intersector_t, local_t, links_t, array_t, 4> {
    public:
    using base_type = mask_base<intersector_t, local_t, links_t, array_t, 4>;
    using base_type::base_type;
    using mask_tolerance = typename base_type::template array_type<scalar, 2>;
    using mask_values = typename base_type::mask_values;
    using links_type = typename base_type::links_type;
    using local_type = typename base_type::local_type;
    using intersector_type = typename base_type::intersector_type;
    using point2 = __plugin::point2<scalar>;

    static constexpr mask_tolerance within_epsilon = {
        std::numeric_limits<scalar>::epsilon(),
        std::numeric_limits<scalar>::epsilon()};

    /* Default constructor */
    trapezoid2()
        : base_type({std::numeric_limits<scalar>::infinity(),
                     std::numeric_limits<scalar>::infinity(),
                     std::numeric_limits<scalar>::infinity(),
                     std::numeric_limits<scalar>::infinity()},
                    {}) {}

    /** Construction from boundary values
     *
     * @param half_length_0 first half length in loc0
     * @param half_length_1 second half length in loc0
     * @param half_length_2 half length in loc1
     */
    DETRAY_HOST_DEVICE
    trapezoid2(scalar half_length_0, scalar half_length_1, scalar half_length_2,
               links_type links)
        : base_type({half_length_0, half_length_1, half_length_2,
                     static_cast<scalar>(1. / (2. * half_length_2))},
                    links) {}

    /** Assignment operator from an array, convenience function
     *
     * @param rhs is the right hand side object
     **/
    DETRAY_HOST_DEVICE
    trapezoid2<intersector_t, local_type, links_type, array_t> &operator=(
        const array_t<scalar, 3> &rhs) {
        this->_values = rhs;
        return (*this);
    }

    /** Mask operation
     *
     * @tparam inside_local_t is the type of the local frame (ignored here)
     *
     * @param p the point to be checked
     * @param t us the tolerance tuple (l0,l1)
     *
     * @return an intersection status e_inside / e_outside
     **/
    template <typename inside_local_t>
    DETRAY_HOST_DEVICE intersection::status is_inside(
        const point2 &p, const mask_tolerance t = within_epsilon) const {
        scalar rel_y = (this->_values[2] + p[1]) * this->_values[3];
        return (std::abs(p[0]) <=
                    this->_values[0] +
                        rel_y * (this->_values[1] - this->_values[0]) + t[0] and
                std::abs(p[1]) <= this->_values[2] + t[1])
                   ? intersection::status::e_inside
                   : intersection::status::e_outside;
    }

    /** Transform to a string for output debugging */
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "trapezoid2";
        for (const auto &v : this->_values) {
            ss << ", " << v;
        }
        return ss.str();
    }
};

}  // namespace detray
