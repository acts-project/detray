/** Detray library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s).
#include <cmath>

namespace detray {

/// @brief Local frame projection into a 2D cylindrical coordinate frame
template <typename algebra_t>
struct cylindrical2D {

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

    /// This method transform from a point from global cartesian 3D frame to a
    /// local 2D cylindrical point
    DETRAY_HOST_DEVICE
    static inline point3D global_to_local(const transform3D &trf,
                                          const point3D &p,
                                          const vector3D & /*d*/) {
        const auto local3 = trf.point_to_local(p);

        return {getter::perp(local3) * getter::phi(local3), local3[2],
                getter::perp(local3)};
    }

    /// This method transform from a local 2D cylindrical point to a point
    /// global cartesian 3D frame
    DETRAY_HOST_DEVICE
    static inline point3D local_to_global(const transform3D &trf,
                                          const point3D &p) {
        const scalar_t r{p[2]};
        const scalar_t phi{p[0] / r};
        const scalar_t x{r * math_ns::cos(phi)};
        const scalar_t y{r * math_ns::sin(phi)};
        const scalar_t z{p[1]};

        return trf.point_to_global(point3D{x, y, z});
    }

    /// This method transform from a local 2D cylindrical point to a point
    /// global cartesian 3D frame
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3D bound_local_to_global(
        const transform3D &trf, const mask_t &mask, const point2D &p,
        const vector3D & /*dir*/) {

        return local_to_global(trf, {p[0], p[1], mask[mask_t::shape::e_r]});
    }

    /// @returns the normal vector given a bound position @param bound_pos
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3D normal(const transform3D &trf3,
                                                     const point2D &bound_pos,
                                                     const mask_t &mask) {
        const scalar_t phi{bound_pos[0] / mask[mask_t::shape::e_r]};
        const vector3D local_normal{math_ns::cos(phi), math_ns::sin(phi), 0.f};

        // normal vector in local coordinate
        return trf3.rotate(trf3.matrix(), local_normal);
    }

    /// @returns the normal vector given a local position @param loc_pos
    DETRAY_HOST_DEVICE
    static inline vector3D normal(const transform3D &trf3,
                                  const point3D &loc_pos) {
        const scalar_t phi{loc_pos[0] / loc_pos[2]};
        const vector3D local_normal{math_ns::cos(phi), math_ns::sin(phi), 0.f};

        // normal vector in local coordinate
        return trf3.rotate(trf3.matrix(), local_normal);
    }

};  // struct cylindrical2

}  // namespace detray
