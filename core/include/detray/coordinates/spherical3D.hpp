/** Detray library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/boolean.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

template <typename algebra_t>
struct spherical3D {

    /// @name Type definitions for the struct
    /// @{
    using scalar_t = dscalar<algebra_t>;
    using bool_t = dbool<algebra_t>;
    using point2D = dpoint2D<algebra_t>;
    using point3D = dpoint3D<algebra_t>;
    using vector3D = dvector3D<algebra_t>;
    using transform3D = dtransform3D<algebra_t>;

    // Local point type on a 2D spherical surface
    using loc_point = point2D;

    /// @brief project to a 2D spherical local frame
    DETRAY_HOST_DEVICE
    constexpr loc_point to_2D_local(const point3D &p) const {
        return {p[1], p[2]};
    }

    /// @brief transform a point from 3D cartesian to 3D sherical frame
    DETRAY_HOST_DEVICE
    point3D global_to_local(const transform3D &trf, const point3D &p,
                            const vector3D & /*d*/) const {
        const auto local3 = trf.point_to_local(p);
        return {getter::norm(local3), getter::phi(local3),
                getter::theta(local3)};
    }

    /// @brief transform a point from 3D sherical to 3D cartesian frame
    DETRAY_HOST_DEVICE
    point3D local_to_global(const transform3D &trf,
                            const point3D &loc_p) const {
        const auto sin_theta{math_ns::sin(loc_p[2])};

        const point3D glob_p{sin_theta * math_ns::cos(loc_p[1]),
                             sin_theta * math_ns::sin(loc_p[1]),
                             math_ns::cos(loc_p[2])};

        return trf.point_to_global(loc_p[0] * glob_p);
    }

    /// @returns the normal vector given a bound position @param bound_pos
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline vector3D normal(const transform3D &,
                                              const point2D &bound_pos,
                                              const mask_t &) const {
        const auto sin_theta{math_ns::sin(bound_pos[1])};

        return {sin_theta * math_ns::cos(bound_pos[0]),
                sin_theta * math_ns::sin(bound_pos[0]),
                math_ns::cos(bound_pos[1])};
    }

    /// @returns the normal vector given a local position @param loc_pos
    DETRAY_HOST_DEVICE inline vector3D normal(const transform3D &,
                                              const point3D &loc_pos) const {
        const auto sin_theta{math_ns::sin(loc_pos[2])};

        return {sin_theta * math_ns::cos(loc_pos[1]),
                sin_theta * math_ns::sin(loc_pos[1]), math_ns::cos(loc_pos[2])};
    }
};

}  // namespace detray
