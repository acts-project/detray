/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <climits>
#include <cmath>
#include <optional>
#include <sstream>
#include <string>

#include "core/intersection.hpp"
#include "masks/mask_identifier.hpp"
#include "tools/cylinder_intersector.hpp"
#include "definitions/detray_qualifiers.hpp"

namespace detray {
/** This is a simple mask for a full cylinder
 *
 * @tparam kRadialCheck is a boolean to steer wheter the radius compatibility
 *needs to be checked
 * @tparam intersector_type is a struct used for intersecting this cylinder
 * @tparam links_type is an object where the mask can link to
 * @tparam kMaskContext is a unique mask identifier in a certain context
 *
 * It is defined by r and the half length.
 *
 * @note  While the mask_context can change depending on the typed container
 * structure the mask_identifier is a const expression that determines the
 * mask type once for all.
 *
 **/
template <bool kRadialCheck = true,
          typename intersector_type = detray::cylinder_intersector,
          typename local_type = __plugin::cylindrical2,
          typename links_type = bool, unsigned int kMaskContext = e_cylinder3,
          template <typename, unsigned int> class array_type = darray>
struct cylinder3 {

    using mask_tolerance = array_type<scalar, 2>;

    // This masks checks on: radius, -z, +z
    using mask_values = array_type<scalar, 3>;

    using mask_links_type = links_type;

    mask_values _values = {std::numeric_limits<scalar>::infinity(),
                           -std::numeric_limits<scalar>::infinity(),
                           std::numeric_limits<scalar>::infinity()};

    links_type _links;

    static constexpr unsigned int mask_context = kMaskContext;

    static constexpr unsigned int mask_identifier = e_cylinder3;

    static constexpr mask_tolerance within_epsilon = {
        std::numeric_limits<scalar>::epsilon(),
        std::numeric_limits<scalar>::epsilon()};

    /** Construction from boundary values
     *
     * @param r radius
     * @param half_length_1 half_length
     * @param half_length_2 half length
     */
    DETRAY_HOST_DEVICE
    cylinder3(scalar r, scalar half_length_1, scalar half_length_2 = 0.)
        : _values{r, half_length_1, half_length_2} {}

    /** Assignment operator from an array, convenience function
     *
     * @param rhs is the right hand side object
     **/
    DETRAY_HOST_DEVICE
    cylinder3<kRadialCheck, intersector_type, local_type, links_type,
              kMaskContext>
        &operator=(const array_type<scalar, 3> &rhs) {
        _values = rhs;
        return (*this);
    }

    /** Mask operation
     *
     * @tparam inside_local_type::point3 is the deduced type of the point to be
     *checked w.r.t. to the mask bounds, it's assumed to be within the cylinder
     *3D frame
     *
     * @param p the point to be checked
     * @param t is the tolerance tuple in (radius, z)
     *
     * @return an intersection status e_inside / e_outside
     **/
    template <typename inside_local_type>
    DETRAY_HOST_DEVICE
    intersection_status is_inside(
        const point3 &p, const mask_tolerance t = within_epsilon) const {
        if constexpr (kRadialCheck) {
            scalar r = getter::perp(p);
            if (std::abs(r - _values[0]) >=
                t[0] + 5 * std::numeric_limits<scalar>::epsilon()) {
                return e_missed;
            }
        }
        return (_values[1] - t[1] <= p[2] and p[2] <= _values[2] + t[1])
                   ? e_inside
                   : e_outside;
    }

    /** Equality operator from an array, convenience function
     *
     * @param rhs is the rectangle to be compared with
     *
     * checks identity within epsilon and @return s a boolean*
     **/
    DETRAY_HOST_DEVICE
    bool operator==(const array_type<scalar, 3> &rhs) {
        return (_values == rhs);
    }

    /** Equality operator
     *
     * @param rhs is the rectangle to be compared with
     *
     * checks identity within epsilon and @return s a boolean*
     **/
    DETRAY_HOST_DEVICE
    bool operator==(const cylinder3<> &rhs) { return operator==(rhs._values); }

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
    intersector_type intersector() const { return intersector_type{}; };

    /** Return the values */
    DETRAY_HOST_DEVICE
    const mask_values &values() const { return _values; }

    /** Return the local frame type */
    DETRAY_HOST_DEVICE
    local_type local() const { return local_type{}; }

    /** Return the volume link - const reference */
    DETRAY_HOST_DEVICE
    const links_type &links() const { return _links; }

    /** Return the volume link - non-const access */
    DETRAY_HOST_DEVICE
    links_type &links() { return _links; }

    /** Transform to a string for output debugging */
    DETRAY_HOST_DEVICE
    std::string to_string() const {
        std::stringstream ss;
        ss << "cylinder3," << kMaskContext;
        for (const auto &v : _values) {
            ss << "," << v;
        }
        return ss.str();
    }
};

}  // namespace detray
