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

#include "detray/core/intersection.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/masks/mask_identifier.hpp"
#include "detray/tools/planar_intersector.hpp"

namespace detray {
/** This is a 2-dimensional mask for the annulus geometry that is
 *  e.g. used for the itk strip endcaps.
 *
 * @tparam intersector_type is a struct used for intersecting this cylinder
 * @tparam local_type is the default local type for this mask
 * @tparam links_type is an object where the mask can link to
 * @tparam kMaskContext is a unique mask identifier in a certain context
 *
 * It is defined by the two radii _values[0] and  _values[1] in the polar
 * coordinate system of and endcap strip module, as well as the two phi
 * boundaries, _values[2] and _values[3], that are only easily definable in
 * the local strips system (relative to the average Phi along the
 * strips r coordinate). Using a conversion between the two coordinate
 * systems, these boundaries can be checked with a tolerance in r (t[0] and
 * t[1]), as well as phi (t[3] and t[3]).
 * Due to the local polar coordinate system of the strips needing a
 * different origin from the discs polar one, three additional conversion
 * parameters are included (_values[4], values[5], _values[6]).
 * Here, the first two are the origin shift in xy, while _values[6] is the
 * average Phi angle mentioned above.
 *
 * @note  While the mask_context can change depending on the typed container
 * structure the mask_identifier is a const expression that determines the
 * mask type once for all.
 *
 **/
template <typename intersector_type = planar_intersector,
          typename mask_local_type = __plugin::polar2<detray::scalar>,
          typename mask_links_type = unsigned int,
          unsigned int kMaskContext = e_annulus2,
          template <typename, unsigned int> class array_type = darray>
struct annulus2 {
    using mask_tolerance = array_type<scalar, 2>;
    using mask_values = array_type<scalar, 7>;
    // Export those types
    using links_type = mask_links_type;
    using local_type = mask_local_type;

    mask_values _values = {0.,
                           std::numeric_limits<scalar>::infinity(),
                           -std::numeric_limits<scalar>::infinity(),
                           std::numeric_limits<scalar>::infinity(),
                           0.,
                           0.,
                           0.};

    links_type _links;

    static constexpr unsigned int mask_context = kMaskContext;

    static constexpr unsigned int mask_identifier = e_annulus2;

    static constexpr mask_tolerance within_epsilon = {
        std::numeric_limits<scalar>::epsilon(),
        std::numeric_limits<scalar>::epsilon()};

    /** Construction from boundary values
     *
     * @param r_low lower r boundary
     * @param r_high upper r boundary
     * @param phi_low lower phi boundary
     * @param phi_high upper phi boundary
     * @param shift_x origin shift loc0
     * @param shift_y origin shift loc1
     * @param avg_phi average phi value
     */
    DETRAY_HOST_DEVICE
    annulus2(scalar r_low, scalar r_high, scalar phi_low, scalar phi_high,
             scalar shift_x = 0., scalar shift_y = 0., scalar avg_phi = 0.)
        : _values{r_low, r_high, phi_low, phi_high, shift_x, shift_y, avg_phi} {
    }

    /** Assignment operator from an array, convenience function
     *
     * @param rhs is the right hand side object
     **/
    DETRAY_HOST_DEVICE
    annulus2<intersector_type, local_type, links_type, kMaskContext> &operator=(
        const array_type<scalar, 7> &rhs) {
        _values = rhs;
        return (*this);
    }

    /** Mask operation
     *
     * @tparam inside_local_type is the local type for checking, needs to be
     * specificed
     *
     * @param p the point to be checked in local polar coord
     * @param t is the tolerance in (r, phi)
     *
     * @return an intersection status e_inside / e_outside
     **/
    template <typename inside_local_type>
    DETRAY_HOST_DEVICE intersection_status
    is_inside(const point2 &p, const mask_tolerance t = within_epsilon) const {
        // The two quantities to check: r^2 in module system, phi in strips
        // system

        // In cartesian coordinates go to modules system by shifting origin
        if constexpr (std::is_same_v<inside_local_type,
                                     __plugin::cartesian2<detray::scalar> >) {
            // Calculate radial coordinate in module system:
            scalar x_mod = p[0] - _values[4];
            scalar y_mod = p[1] - _values[5];
            scalar r_mod2 = x_mod * x_mod + y_mod * y_mod;

            // apply tolerances
            scalar minR_tol = _values[0] - t[0];
            scalar maxR_tol = _values[1] + t[0];

            if (r_mod2 < minR_tol * minR_tol or r_mod2 > maxR_tol * maxR_tol)
                return e_outside;

            scalar phi_strp = getter::phi(p) - _values[6];
            // Check phi boundaries, which are well def. in local frame
            return (phi_strp >= _values[2] - t[1] and
                    phi_strp <= _values[3] + t[1])
                       ? e_inside
                       : e_outside;
        }
        // polar strip coordinates given
        else {
            // For a point p in local polar coordinates, rotate by avr phi
            scalar phi_strp = p[1] - _values[6];

            // Check phi boundaries, which are well def. in local frame
            if (phi_strp < _values[2] - t[1] || phi_strp > _values[3] + t[1])
                return e_outside;

            // Now go to module frame to check r boundaries. Use the origin
            // shift in polar coordinates for that
            point2 shift_xy = {-1 * _values[4], -1 * _values[5]};
            scalar shift_r = getter::perp(shift_xy);
            scalar shift_phi = getter::phi(shift_xy);

            scalar r_mod2 = shift_r * shift_r + p[0] * p[0] +
                            2 * shift_r * p[0] * std::cos(phi_strp - shift_phi);

            // Apply tolerances
            scalar minR_tol = _values[0] - t[0];
            scalar maxR_tol = _values[1] + t[0];

            return (r_mod2 >= minR_tol * minR_tol and
                    r_mod2 <= maxR_tol * maxR_tol)
                       ? e_inside
                       : e_outside;
        }
    }

    /** Equality operator from an array, convenience function
     *
     * @param rhs is the rectangle to be compared with
     *
     **/
    DETRAY_HOST_DEVICE
    bool operator==(const array_type<scalar, 7> &rhs) {
        return (_values == rhs);
    }

    /** Equality operator
     *
     * @param rhs is the rectangle to be compared with
     *
     **/
    DETRAY_HOST_DEVICE
    bool operator==(const annulus2 &rhs) {
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
    intersector_type intersector() const { return intersector_type{}; };

    /** Return the values */
    DETRAY_HOST_DEVICE
    const mask_values &values() const { return _values; }

    /** Return the local frame type */
    DETRAY_HOST_DEVICE
    constexpr local_type local() const { return local_type{}; }

    /** Return the volume link - const reference */
    DETRAY_HOST_DEVICE
    const links_type &links() const { return _links; }

    /** Return the volume link - non-const access */
    DETRAY_HOST_DEVICE
    links_type &links() { return _links; }

    /** Transform to a string for output debugging */
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "annulus2," << kMaskContext;
        for (const auto &v : _values) {
            ss << "," << v;
        }
        return ss.str();
    }
};

}  // namespace detray
