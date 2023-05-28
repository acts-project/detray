/** Detray library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

template <typename transform3_t>
struct spherical3D {

    /// @name Type definitions for the struct
    /// @{

    using transform3_type = transform3_t;
    using scalar_type = typename transform3_t::scalar_type;
    using point2 = typename transform3_t::point2;
    using point3 = typename transform3_t::point3;
    using vector3 = typename transform3_t::vector3;

    // Local point type on a 2D spherical surface
    using loc_point = point2;

    /// @brief project to a 2D spherical local frame
    DETRAY_HOST_DEVICE
    constexpr loc_point to_2D_local(const point3 &p) const {
        return {p[1], p[2]};
    }

    /// @brief transform a point from 3D cartesian to 3D sherical frame
    DETRAY_HOST_DEVICE
    point3 global_to_local(const transform3_t &trf, const point3 &p,
                           const vector3 & /*d*/) const {
        const auto local3 = trf.point_to_local(p);
        return {getter::perp(local3), getter::phi(local3), getter::theta(local3)};
    }

    /// @brief transform a point from 3D sherical to 3D cartesian frame
    DETRAY_HOST_DEVICE 
    point3 local_to_global(const transform3_t &trf, 
                           const point3 &loc_p) const {
        const scalar_type sin_theta{math_ns::sin(loc_p[1])};

        const point3 glob_p{sin_theta * math_ns::cos(loc_p[1]),
                            sin_theta * math_ns::sin(loc_p[1]),
                            math_ns::cos(loc_p[2])};

        return trf.point_to_global(loc_p[0] * glob_p);
    }

    /// This method transform from a local 2D polar point to a point global
    /// cartesian 3D frame
    /*template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 bound_local_to_global(
        const transform3_t &trf, const mask_t &, const point2 &p,
        const vector3 & ) const {

        return this->local_to_global(trf, {p[0], p[1], 0.f});
    }*/

    DETRAY_HOST_DEVICE 
    vector3 normal(const transform3_t &trf3,
                   const point3 & loc_p) const {
        return local_to_global(trf3, vector::normalize(loc_p));
    }
};

}  // namespace detray
