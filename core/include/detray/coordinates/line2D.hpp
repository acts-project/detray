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

namespace detray {

template <typename algebra_t>
struct line2D {

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
    /// local 2D line point
    DETRAY_HOST_DEVICE
    static inline point3D global_to_local(const transform3D &trf,
                                          const point3D &p, const vector3D &d) {

        const auto local3 = trf.point_to_local(p);

        // Line direction
        const vector3D z = getter::vector<3>(trf.matrix(), 0u, 2u);

        // Line center
        const point3D tp = trf.translation() - p;

        // Radial vector
        const vector3D r = vector::cross(z, d);

        // Assign the sign depending on the position w.r.t line
        // Right: -1
        // Left: 1
        const scalar_t sign =
            math_ns::copysign(scalar_t(1.f), vector::dot(r, tp));

        return {sign * getter::perp(local3), local3[2], getter::phi(local3)};
    }

    /// This method transform from a local 2D line point to a point global
    /// cartesian 3D frame
    DETRAY_HOST_DEVICE
    static inline point3D local_to_global(const transform3D &trf,
                                          const point3D &p) {
        const scalar_t R = math_ns::abs(p[0]);
        const point3D local = {R * math_ns::cos(p[2]), R * math_ns::sin(p[2]),
                               p[1]};

        return trf.point_to_global(local);
    }

    /// This method transform from a local 2D line point to a point global
    /// cartesian 3D frame
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3D bound_local_to_global(
        const transform3D &trf, const mask_t & /*mask*/, const point2D &p,
        const vector3D &d) {

        // Line direction
        const vector3D z = getter::vector<3>(trf.matrix(), 0u, 2u);

        // Radial vector
        const vector3D r = vector::cross(z, d);

        // Local Z poisition in global cartesian coordinate
        const point3D locZ_in_global =
            trf.point_to_global(point3D{0.f, 0.f, p[1]});

        return locZ_in_global + p[0] * vector::normalize(r);
    }

    /// @returns the normal vector
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3D normal(const transform3D &trf3,
                                                     const point2D & = {},
                                                     const mask_t & = {}) {
        return getter::vector<3>(trf3.matrix(), 0u, 2u);
    }

    /// @returns the normal vector
    DETRAY_HOST_DEVICE
    static inline vector3D normal(const transform3D &trf3,
                                  const point3D & = {}) {
        return getter::vector<3>(trf3.matrix(), 0u, 2u);
    }
};

}  // namespace detray
