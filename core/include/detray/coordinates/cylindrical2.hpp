/** Detray library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinate_base.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s).
#include <cmath>

namespace detray {

/** Local frame projection into a 2D cylindrical coordinate frame
 */
template <typename transform3_t>
struct cylindrical2 : public coordinate_base<cylindrical2, transform3_t> {

    /// @name Type definitions for the struct
    /// @{

    /// Base type
    using base_type = coordinate_base<cylindrical2, transform3_t>;
    /// Sclar type
    using scalar_type = typename transform3_t::scalar_type;
    /// Point in 2D space
    using point2 = typename transform3_t::point2;
    /// Point in 3D space
    using point3 = typename transform3_t::point3;
    /// Vector in 3D space
    using vector3 = typename transform3_t::vector3;

    /// @}

    /** This method transform from a point from 3D cartesian frame to a 2D
     * cylindrical point */
    DETRAY_HOST_DEVICE
    inline point2 operator()(const point3 &p) const {

        return {getter::perp(p) * getter::phi(p), p[2]};
    }

    /** This method transform from a point from global cartesian 3D frame to a
     * local 2D cylindrical point */
    DETRAY_HOST_DEVICE
    inline point2 global_to_local(const transform3_t &trf, const point3 &p,
                                  const vector3 & /*d*/) const {
        const auto local3 = trf.point_to_local(p);
        return this->operator()(local3);
    }

    /** This method transform from a local 2D cylindrical point to a point
     * global cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 local_to_global(
        const transform3_t &trf, const mask_t &mask, const point2 &p,
        const vector3 & /*d*/) const {
        const scalar_type r = mask[0];
        const scalar_type phi = p[0] / r;
        const scalar_type x = r * std::cos(phi);
        const scalar_type y = r * std::sin(phi);
        const scalar_type z = p[1];

        return trf.point_to_global(point3{x, y, z});
    }

};  // struct cylindrical2

}  // namespace detray