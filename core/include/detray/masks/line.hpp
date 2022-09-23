/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinates.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/line_intersector.hpp"
#include "detray/masks/mask_base.hpp"

// System include(s).
#include <climits>
#include <cmath>
#include <type_traits>

namespace detray {

/** This is a simple mask for a line defined with a line length and its scope
 *
 **/
template <typename transform3_t = __plugin::transform3<scalar>,
          template <class> typename intersector_t = line_intersector,
          template <class> typename local_t = cartesian2,
          typename links_t = dindex, bool kSquareScope = false,
          template <typename, std::size_t> class array_t = darray>
class line final : public mask_base<transform3_t, intersector_t, local_t,
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
    using point3 = typename transform3_t::point3;

    static constexpr bool square_scope = kSquareScope;

    /* Default constructor */
    line() : base_type({0, std::numeric_limits<scalar>::infinity()}, {}) {}

    /** Construction from boundary values
     *
     * @param half_cell_size is the half length of transverse cell of line
     * detector
     * @param half_length is the half length of line
     */
    DETRAY_HOST_DEVICE
    line(scalar half_cell_size, scalar half_length, links_type links)
        : base_type({half_cell_size, half_length}, links) {}

    /** Mask operation
     *
     * @tparam inside_local_t is the local type for inside checking
     *
     * @param p is the intersection point in local cartesian coordinate
     * @param t is the tolerance in r
     *
     * @return an intersection status e_inside / e_outside
     **/
    DETRAY_HOST_DEVICE intersection::status is_inside(
        const point3 &p,
        const scalar t = std::numeric_limits<scalar>::epsilon()) const {

        // For square cross section, we check if (1) the x and y of transverse
        // position is less than the half cell size and (2) the distance to the
        // point of closest approach on line from the line center is less than
        // the half line length
        if constexpr (square_scope) {
            return std::abs(p[0]) <= this->_values[0] + t &&
                           std::abs(p[1]) <= this->_values[0] + t &&
                           std::abs(p[2]) <= this->_values[1] + t
                       ? intersection::status::e_inside
                       : intersection::status::e_outside;

            // For circular cross section, we check if (1) the radial distance
            // is within the scope and (2) the distance to the point of closest
            // approach on line from the line center is less than the half line
            // length
        } else {
            return (getter::perp(p) <= this->_values[0] + t &&
                    std::abs(p[2]) <= this->_values[1] + t)

                       ? intersection::status::e_inside
                       : intersection::status::e_outside;
        }
    }

    /** Transform to a string for output debugging */
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "line";
        for (const auto &v : this->_values) {
            ss << ", " << v;
        }
        return ss.str();
    }
};

}  // namespace detray