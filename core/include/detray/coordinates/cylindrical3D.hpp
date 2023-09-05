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

/// @brief Frame projection into a cartesian coordinate frame
template <typename algebra_t>
struct cylindrical3D {

    /// @name Linear algebra definitions
    /// @{
    using algebra = algebra_t;
    using scalar_t = dscalar<algebra>;
    using point3D = dpoint3D<algebra>;
    using vector3D = dvector3D<algebra>;
    using transform3D = dtransform3D<algebra>;

    // Local point type on a 2D spherical surface
    using loc_point = point3D;
    /// @}

    /// This method transform from a point from global cartesian 3D frame to a
    /// local 3D cylindrical point
    DETRAY_HOST_DEVICE
    static inline point3D global_to_local(const transform3D &trf,
                                          const point3D &p,
                                          const vector3D & /*d*/) {
        const auto local3 = trf.point_to_local(p);
        return {getter::perp(local3), getter::phi(local3), local3[2]};
    }

    /// This method transform from a local 3D cylindrical point to a point
    /// global cartesian 3D frame
    DETRAY_HOST_DEVICE
    static inline point3D local_to_global(const transform3D &trf,
                                          const point3D &p) {
        const scalar_t x{p[0] * math_ns::cos(p[1])};
        const scalar_t y{p[0] * math_ns::sin(p[1])};

        return trf.point_to_global(point3D{x, y, p[2]});
    }

};  // struct cylindrical3

}  // namespace detray
