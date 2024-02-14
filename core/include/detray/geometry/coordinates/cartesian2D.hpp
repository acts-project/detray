/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {

/// Projection into a 2D cartesian coordinate frame
template <typename algebra_t>
struct cartesian2D {

    using transform3_type = algebra_t;
    using scalar_type = typename algebra_t::scalar_type;
    using point2 = typename algebra_t::point2;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    /// Local point type in 2D cartesian coordinates
    using loc_point = point2;

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 3D cartesian point
    DETRAY_HOST_DEVICE
    static inline point3 global_to_local_3D(const transform3_type &trf,
                                            const point3 &p,
                                            const vector3 & /*dir*/) {
        return trf.point_to_local(p);
    }

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 2D cartesian point
    DETRAY_HOST_DEVICE
    static inline loc_point global_to_local(const transform3_type &trf,
                                            const point3 &p,
                                            const vector3 & /*dir*/) {
        auto loc_p = trf.point_to_local(p);
        return {loc_p[0], loc_p[1]};
    }

    /// This method transforms from a local 3D cartesian point to a point in
    /// the global cartesian 3D frame
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_type &trf, const point3 &p) {
        return trf.point_to_global(p);
    }

    /// This method transforms from a local 2D cartesian point to a point in
    /// the global cartesian 3D frame
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_type &trf, const mask_t & /*mask*/, const loc_point &p,
        const vector3 & /*dir*/) {
        return trf.point_to_global(point3{p[0], p[1], 0.f});
    }

    /// @returns the normal vector in global coordinates
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_type &trf,
                                                    const point2 & = {},
                                                    const mask_t & = {}) {
        return trf.z();
    }

    /// @returns the normal vector in global coordinates
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_type &trf,
                                                    const point3 & = {}) {
        return trf.z();
    }

};  // struct cartesian2D

}  // namespace detray
