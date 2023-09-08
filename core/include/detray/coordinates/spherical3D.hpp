/** Detray library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

template <typename algebra_t>
struct spherical3D {

    /// @name Linear algebra definitions
    /// @{
    using algebra = algebra_t;
    using scalar_t = dscalar<algebra>;
    using point2D = dpoint2D<algebra>;
    using point3D = dpoint3D<algebra>;
    using vector3D = dvector3D<algebra>;
    using transform3D = dtransform3D<algebra>;

    // Local point type on a 2D spherical surface
    using loc_point = point2D;
    /// @}

    /// @brief transform a point from 3D cartesian to 2D sherical frame
    DETRAY_HOST_DEVICE
    static inline point2D global_to_bound(const transform3D &trf,
                                          const point3D &p,
                                          const vector3D & /*d*/) {
        const auto local3 = trf.point_to_local(p);
        return {getter::phi(local3), getter::theta(local3)};
    }

    /// @brief transform a point from 3D cartesian to 3D sherical frame
    DETRAY_HOST_DEVICE
    static inline point3D global_to_local(const transform3D &trf,
                                          const point3D &p,
                                          const vector3D & /*d*/) {
        const auto local3 = trf.point_to_local(p);
        return {getter::norm(local3), getter::phi(local3),
                getter::theta(local3)};
    }

    /// @returns the global cartesian point from a bound local (2D) point
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3D bound_local_to_global(
        const transform3D &trf3, const mask_t &mask, const point2D &p,
        const vector3D & /*d*/) {

        return local_to_global(trf3, {mask[mask_t::shape::e_r], p[0], p[1]});
    }

    /// @brief transform a point from 3D sherical to 3D cartesian frame
    DETRAY_HOST_DEVICE
    static inline point3D local_to_global(const transform3D &trf,
                                          const point3D &loc_p) {
        const auto sin_theta{math_ns::sin(loc_p[2])};

        const point3D glob_p{sin_theta * math_ns::cos(loc_p[1]),
                             sin_theta * math_ns::sin(loc_p[1]),
                             math_ns::cos(loc_p[2])};

        return trf.point_to_global(loc_p[0] * glob_p);
    }

    /// @returns the normal vector given a bound position @param bound_pos
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3D normal(const transform3D &,
                                                     const point2D &bound_pos,
                                                     const mask_t &) {
        const auto sin_theta{math_ns::sin(bound_pos[1])};

        return {sin_theta * math_ns::cos(bound_pos[0]),
                sin_theta * math_ns::sin(bound_pos[0]),
                math_ns::cos(bound_pos[1])};
    }

    /// @returns the normal vector given a local position @param loc_pos
    DETRAY_HOST_DEVICE
    static inline vector3D normal(const transform3D &, const point3D &loc_pos) {
        const auto sin_theta{math_ns::sin(loc_pos[2])};

        return {sin_theta * math_ns::cos(loc_pos[1]),
                sin_theta * math_ns::sin(loc_pos[1]), math_ns::cos(loc_pos[2])};
    }
};

}  // namespace detray
