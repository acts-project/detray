/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical2.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/masks/mask_base.hpp"

// System include(s)
#include <climits>
#include <cmath>
#include <sstream>
#include <string>

namespace detray {
/** This is a simple mask for a full cylinder
 *
 * @tparam kRadialCheck is a boolean to steer wheter the radius compatibility
 *needs to be checked
 * @tparam intersector_t is a struct used for intersecting this cylinder
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
template <typename transform3_t = __plugin::transform3<scalar>,
          template <class> typename intersector_t = cylinder_intersector,
          template <class> typename local_t = cylindrical2,
          typename links_t = dindex, bool kRadialCheck = false,
          template <typename, std::size_t> class array_t = darray>
class cylinder3 final : public mask_base<transform3_t, intersector_t, local_t,
                                         links_t, array_t, 3> {
    public:
    using base_type =
        mask_base<transform3_t, intersector_t, local_t, links_t, array_t, 3>;
    using base_type::base_type;
    using mask_values = typename base_type::mask_values;
    using links_type = typename base_type::links_type;
    using local_type = typename base_type::local_type;
    using intersector_type = typename base_type::intersector_type;
    using point3 = typename transform3_t::point3;
    using scalar_type = typename transform3_t::scalar_type;

    scalar_type m_r;

    /* Default constructor */
    cylinder3()
        : base_type({std::numeric_limits<scalar>::infinity(),
                     -std::numeric_limits<scalar>::infinity(),
                     std::numeric_limits<scalar>::infinity()},
                    {}) {}

    /** Construction from boundary values
     *
     * @param r radius
     * @param half_length_1 half_length
     * @param half_length_2 half length
     */
    DETRAY_HOST_DEVICE
    cylinder3(scalar r, scalar half_length_1, scalar half_length_2,
              links_type links)
        : base_type({r, half_length_1, half_length_2}, links), m_r(r) {}

    /// Get radius
    DETRAY_HOST_DEVICE
    scalar_type radius() const { return m_r; }

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
    template <typename inside_local_t, bool is_rad_check = kRadialCheck>
    DETRAY_HOST_DEVICE intersection::status is_inside(
        const point3 &p,
        const scalar t = std::numeric_limits<scalar>::epsilon()) const {
        if constexpr (is_rad_check) {
            scalar r = getter::perp(p);
            if (std::abs(r - this->_values[0]) >=
                t + 5 * std::numeric_limits<scalar>::epsilon()) {
                return intersection::status::e_missed;
            }
        }
        return (this->_values[1] - t <= p[2] and p[2] <= this->_values[2] + t)
                   ? intersection::status::e_inside
                   : intersection::status::e_outside;
    }

    /** Transform to a string for output debugging */
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream ss;
        ss << "cylinder3";
        for (const auto &v : this->_values) {
            ss << ", " << v;
        }
        return ss.str();
    }
};

}  // namespace detray
