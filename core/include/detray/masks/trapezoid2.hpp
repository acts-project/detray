/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <climits>
#include <cmath>
#include <sstream>
#include <string>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/planar_intersector.hpp"

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
template <typename intersector_t = planar_intersector,
          typename mask_local_t = __plugin::cartesian2<detray::scalar>,
          typename mask_links_t = dindex,
          template <typename, std::size_t> class array_t = darray>
struct trapezoid2 {
    using mask_tolerance = array_t<scalar, 2>;
    using mask_values = array_t<scalar, 4>;
    using links_type = mask_links_t;
    using local_type = mask_local_t;

    mask_values _values = {std::numeric_limits<scalar>::infinity(),
                           std::numeric_limits<scalar>::infinity(),
                           std::numeric_limits<scalar>::infinity(),
                           std::numeric_limits<scalar>::infinity()};

    links_type _links;

    static constexpr mask_tolerance within_epsilon = {
        std::numeric_limits<scalar>::epsilon(),
        std::numeric_limits<scalar>::epsilon()};

    /* Default constructor */
    trapezoid2() = default;

    /** Construction from boundary values
     *
     * @param half_length_0 first half length in loc0
     * @param half_length_1 second half length in loc0
     * @param half_length_2 half length in loc1
     */
    DETRAY_HOST_DEVICE
    trapezoid2(scalar half_length_0, scalar half_length_1, scalar half_length_2,
               links_type links)
        : _values{half_length_0, half_length_1, half_length_2,
                  static_cast<scalar>(1. / (2. * half_length_2))},
          _links(links) {}

    /** Assignment operator from an array, convenience function
     *
     * @param rhs is the right hand side object
     **/
    DETRAY_HOST_DEVICE
    trapezoid2<intersector_t, local_type, links_type> &operator=(
        const array_t<scalar, 3> &rhs) {
        _values = rhs;
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
        scalar rel_y = (_values[2] + p[1]) * _values[3];
        return (std::abs(p[0]) <=
                    _values[0] + rel_y * (_values[1] - _values[0]) + t[0] and
                std::abs(p[1]) <= _values[2] + t[1])
                   ? intersection::status::e_inside
                   : intersection::status::e_outside;
    }

    /** Equality operator from an array, convenience function
     *
     * @param rhs is the rectangle to be compared with
     *
     * checks identity within epsilon and @return s a boolean*
     **/
    DETRAY_HOST_DEVICE
    bool operator==(const array_t<scalar, 3> &rhs) { return (_values == rhs); }

    /** Equality operator
     *
     * @param rhs is the rectangle to be compared with
     *
     * checks identity within epsilon and @return s a boolean*
     **/
    DETRAY_HOST_DEVICE
    bool operator==(const trapezoid2 &rhs) {
        return (_values == rhs._values && _links == rhs._links);
    }

    /** Access operator - non-const
     * @return the reference to the member variable
     */
    DETRAY_HOST_DEVICE
    scalar &operator[](unsigned int value_index) {
        return _values[value_index];
    }

    /** Access operator - non-const
     * @return a copy of the member variable
     */
    DETRAY_HOST_DEVICE
    scalar operator[](unsigned int value_index) const {
        return _values[value_index];
    }

    /** Return an associated intersector type */
    DETRAY_HOST_DEVICE
    intersector_t intersector() const { return intersector_t{}; };

    /** Return the values */
    DETRAY_HOST_DEVICE
    const mask_values &values() const { return _values; }

    /** Return the local frame type */
    DETRAY_HOST_DEVICE
    constexpr local_type local() const { return local_type{}; }

    /** @return the links - const reference */
    DETRAY_HOST_DEVICE
    const links_type &links() const { return _links; }

    /** @return the links - non-const access */
    DETRAY_HOST_DEVICE
    links_type &links() { return _links; }

    /** @return the volume link - const reference */
    DETRAY_HOST_DEVICE
    dindex volume_link() const { return detail::get<0>(_links); }

    /** @return the volume link - non-const access */
    DETRAY_HOST_DEVICE
    dindex volume_link() { return detail::get<0>(_links); }

    /** @return the surface finder link - const reference */
    DETRAY_HOST_DEVICE
    dindex finder_link() const { return detail::get<1>(_links); }

    /** @return the surface finder link - non-const access */
    DETRAY_HOST_DEVICE
    dindex finder_link() { return detail::get<1>(_links); }

    /** Transform to a string for output debugging */
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "trapezoid2";
        for (const auto &v : _values) {
            ss << ", " << v;
        }
        return ss.str();
    }
};

}  // namespace detray
