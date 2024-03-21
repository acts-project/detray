/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {

/// Projection into a 3D cylindrical coordinate frame
template <typename algebra_t>
struct cylindrical3D {

    using transform3_type = algebra_t;
    using scalar_type = typename algebra_t::scalar_type;
    using point2 = typename algebra_t::point2;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    /// Local point type in 3D cylindrical coordinates
    using loc_point = point3;

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 3D cylindrical point
    DETRAY_HOST_DEVICE
    static inline point3 global_to_local_3D(const transform3_type &trf,
                                            const point3 &p,
                                            const vector3 &dir) {
        return cylindrical3D<algebra_t>::global_to_local(trf, p, dir);
    }

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 3D cylindrical point
    DETRAY_HOST_DEVICE
    static inline loc_point global_to_local(const transform3_type &trf,
                                            const point3 &p,
                                            const vector3 & /*dir*/) {
        const auto local3 = trf.point_to_local(p);
        return {getter::perp(local3), getter::phi(local3), local3[2]};
    }

    /// This method transforms from a local 3D cylindrical point to a point in
    /// the global cartesian 3D frame
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_type &trf, const point3 &p) {
        const scalar_type x{p[0] * math::cos(p[1])};
        const scalar_type y{p[0] * math::sin(p[1])};

        return trf.point_to_global(point3{x, y, p[2]});
    }

    /// This method transforms from a local 3D cylindrical point to a point in
    /// the global cartesian 3D frame
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_type &trf, const mask_t & /*mask*/, const loc_point &p,
        const vector3 & /*dir*/) {
        return cylindrical3D<algebra_t>::global_to_local(trf, p);
    }

};  // struct cylindrical3D

}  // namespace detray
