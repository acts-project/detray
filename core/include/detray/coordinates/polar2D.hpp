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

template <typename algebra_t>
struct polar2D {

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
    /// local 2D polar point
    DETRAY_HOST_DEVICE
    static inline point3D global_to_local(const transform3D &trf,
                                          const point3D &p,
                                          const vector3D & /*d*/) {
        const auto local3 = trf.point_to_local(p);
        return {getter::perp(local3), getter::phi(local3), local3[2]};
    }

    /// This method transform from a local 2D polar point to a point global
    /// cartesian 3D frame
    DETRAY_HOST_DEVICE
    static inline point3D local_to_global(const transform3D &trf,
                                          const point3D &p) {
        const scalar_t x = p[0] * math_ns::cos(p[1]);
        const scalar_t y = p[0] * math_ns::sin(p[1]);

        return trf.point_to_global(point3D{x, y, p[2]});
    }

    /// This method transform from a local 2D polar point to a point global
    /// cartesian 3D frame
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3D bound_local_to_global(
        const transform3D &trf, const mask_t & /*mask*/, const point2D &p,
        const vector3D & /*d*/) {

        return local_to_global(trf, {p[0], p[1], 0.f});
    }

    /// @returns the normal vector
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3D normal(const transform3D &trf3,
                                                     const point2D & = {},
                                                     const mask_t & = {}) {
        return trf3.z();
    }

    /// @returns the normal vector
    DETRAY_HOST_DEVICE
    static inline vector3D normal(const transform3D &trf3,
                                  const point3D & = {}) {
        return trf3.z();
    }
};

}  // namespace detray
