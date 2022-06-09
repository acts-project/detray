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
#include "detray/masks/mask_base.hpp"

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
template <typename intersector_t = planar_intersector,
          typename local_t = __plugin::cartesian2<detray::scalar>,
          typename links_t = dindex,
          template <typename, std::size_t> class array_t = darray>
class ring2 final : public mask_base<intersector_t, local_t, links_t, array_t> {
    public:
    using base_type = mask_base<intersector_t, local_t, links_t, array_t>;
    using mask_tolerance = scalar;
    using mask_values = typename base_type::template array_type<scalar, 2>;
    using links_type = typename base_type::links_type;
    using local_type = typename base_type::local_type;
    using intersector_type = typename base_type::intersector_type;
    using point2 = __plugin::point2<scalar>;

    static constexpr mask_tolerance within_epsilon =
        std::numeric_limits<scalar>::epsilon();

    /* Default constructor */
    ring2() = default;

    /** Construction from boundary values
     *
     * @param r_low lower radial bound
     * @param r_high upper radial bound
     */
    DETRAY_HOST_DEVICE
    ring2(scalar r_low, scalar r_high, links_type links)
        : _values{r_low, r_high} {
        this->_links = links;
    }

    /** Assignment operator from an array, convenience function
     *
     * @param rhs is the right hand side object
     **/
    DETRAY_HOST_DEVICE
    ring2<intersector_t, local_type, links_type> &operator=(
        const array_t<scalar, 2> &rhs) {
        _values = rhs;
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
        if constexpr (std::is_same_v<inside_local_t,
                                     __plugin::cartesian2<detray::scalar>>) {
            scalar r = getter::perp(p);
            return (r + t >= _values[0] and r <= _values[1] + t)
                       ? intersection::status::e_inside
                       : intersection::status::e_outside;
        }

        return (p[0] + t >= _values[0] and p[0] <= _values[1] + t)
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
    bool operator==(const array_t<scalar, 2> &rhs) { return (_values == rhs); }

    /** Equality operator
     *
     * @param rhs is the rectangle to be compared with
     *
     * checks identity within epsilon and @return s a boolean*
     **/
    DETRAY_HOST_DEVICE
    bool operator==(const ring2 &rhs) {
        return (_values == rhs._values && this->_links == rhs._links);
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

    /** Return the values */
    DETRAY_HOST_DEVICE
    const mask_values &values() const { return _values; }

    /** Transform to a string for output debugging */
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "ring2";
        for (const auto &v : _values) {
            ss << ", " << v;
        }
        return ss.str();
    }

    private:
    mask_values _values = {0., std::numeric_limits<scalar>::infinity()};
};

}  // namespace detray
