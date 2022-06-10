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
#include "detray/intersection/planar_intersector.hpp"
#include "detray/masks/mask_base.hpp"

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

template <typename intersector_t = planar_intersector,
          typename local_t = __plugin::polar2<detray::scalar>,
          typename links_t = dindex,
          template <typename, std::size_t> class array_t = darray>
class annulus2 final
    : public mask_base<intersector_t, local_t, links_t, array_t, 7> {
    public:
    using base_type = mask_base<intersector_t, local_t, links_t, array_t, 7>;
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
    annulus2() = default;

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
             scalar shift_x, scalar shift_y, scalar avg_phi, links_type links)
        : base_type(
              {r_low, r_high, phi_low, phi_high, shift_x, shift_y, avg_phi},
              links) {}

    /** Assignment operator from an array, convenience function
     *
     * @param rhs is the right hand side object
     **/
    DETRAY_HOST_DEVICE
    annulus2<intersector_t, local_type, links_type, array_t> &operator=(
        const mask_values &rhs) {
        this->_values = rhs;
        return (*this);
    }

    /** Mask operation
     *
     * @tparam inside_local_t is the local type for checking, needs to be
     * specificed
     *
     * @param p the point to be checked in local polar coord
     * @param t is the tolerance in (r, phi)
     *
     * @return an intersection status e_inside / e_outside
     **/
    template <typename inside_local_t>
    DETRAY_HOST_DEVICE intersection::status is_inside(
        const point2 &p, const mask_tolerance t = within_epsilon) const {
        // The two quantities to check: r^2 in module system, phi in strips
        // system

        // In cartesian coordinates go to modules system by shifting origin
        if constexpr (std::is_same_v<inside_local_t,
                                     __plugin::cartesian2<detray::scalar>>) {
            // Calculate radial coordinate in module system:
            scalar x_mod = p[0] - this->_values[4];
            scalar y_mod = p[1] - this->_values[5];
            scalar r_mod2 = x_mod * x_mod + y_mod * y_mod;

            // apply tolerances
            scalar minR_tol = this->_values[0] - t[0];
            scalar maxR_tol = this->_values[1] + t[0];

            if (r_mod2 < minR_tol * minR_tol or r_mod2 > maxR_tol * maxR_tol)
                return intersection::status::e_outside;

            scalar phi_strp = getter::phi(p) - this->_values[6];
            // Check phi boundaries, which are well def. in local frame
            return (phi_strp >= this->_values[2] - t[1] and
                    phi_strp <= this->_values[3] + t[1])
                       ? intersection::status::e_inside
                       : intersection::status::e_outside;
        }
        // polar strip coordinates given
        else {
            // For a point p in local polar coordinates, rotate by avr phi
            scalar phi_strp = p[1] - this->_values[6];

            // Check phi boundaries, which are well def. in local frame
            if (phi_strp < this->_values[2] - t[1] ||
                phi_strp > this->_values[3] + t[1])
                return intersection::status::e_outside;

            // Now go to module frame to check r boundaries. Use the origin
            // shift in polar coordinates for that
            point2 shift_xy = {-1 * this->_values[4], -1 * this->_values[5]};
            scalar shift_r = getter::perp(shift_xy);
            scalar shift_phi = getter::phi(shift_xy);

            scalar r_mod2 = shift_r * shift_r + p[0] * p[0] +
                            2 * shift_r * p[0] * std::cos(phi_strp - shift_phi);

            // Apply tolerances
            scalar minR_tol = this->_values[0] - t[0];
            scalar maxR_tol = this->_values[1] + t[0];

            return (r_mod2 >= minR_tol * minR_tol and
                    r_mod2 <= maxR_tol * maxR_tol)
                       ? intersection::status::e_inside
                       : intersection::status::e_outside;
        }
    }

    /// Transform to a string for output debugging
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "annulus2";
        for (const auto &v : this->_values) {
            ss << ", " << v;
        }
        return ss.str();
    }
};

}  // namespace detray
