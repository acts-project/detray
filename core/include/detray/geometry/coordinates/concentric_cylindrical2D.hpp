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

/// Projection into a 2D concentric cylindrical frame
/// (No rotation in coordinate transformation)
template <typename algebra_t>
struct concentric_cylindrical2D {

    using transform3_type = algebra_t;
    using scalar_type = typename algebra_t::scalar_type;
    using point2 = typename algebra_t::point2;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    /// Local point type in 2D cylindrical coordinates
    using loc_point = point2;

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 2D cylindrical point
    DETRAY_HOST_DEVICE
    static inline point3 global_to_local_3D(const transform3_type &trf,
                                            const point3 &p,
                                            const vector3 & /*dir*/) {
        const point3 local3 = p - trf.translation();
        const scalar_type r{getter::perp(local3)};

        return {r * getter::phi(local3), local3[2], r};
    }

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 2D cylindrical point
    DETRAY_HOST_DEVICE
    static inline loc_point global_to_local(const transform3_type &trf,
                                            const point3 &p,
                                            const vector3 & /*dir*/) {
        const point3 local3 = p - trf.translation();

        return {getter::phi(local3), local3[2]};
    }

    /// This method transforms from a local 3D cylindrical point to a point in
    /// the global cartesian 3D frame
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_type &trf, const point3 &p) {

        const scalar_type r{p[2]};
        const scalar_type phi{p[0] / r};
        const scalar_type x{r * math::cos(phi)};
        const scalar_type y{r * math::sin(phi)};
        const scalar_type z{p[1]};

        return point3{x, y, z} + trf.translation();
    }

    /// This method transforms from a local 2D cylindrical point to a point in
    /// the global cartesian 3D frame
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_type &trf, const mask_t &mask, const loc_point &p,
        const vector3 & /*dir*/) {

        const scalar_type r{mask[mask_t::shape::e_r]};
        const scalar_type x{r * math::cos(p[0])};
        const scalar_type y{r * math::sin(p[0])};
        const scalar_type z{p[1]};

        return point3{x, y, z} + trf.translation();
    }

    /// @returns the normal vector in global coordinates
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_type &,
                                                    const point2 &p,
                                                    const mask_t & /*mask*/) {

        // normal vector in global coordinates (concentric cylinders have no
        // rotation)
        return {math::cos(p[0]), math::sin(p[0]), 0.f};
    }

    /// @returns the normal vector given a local position @param p
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_type &,
                                                    const point3 &p) {
        const scalar_type phi{p[0] / p[2]};
        // normal vector in global coordinates (concentric cylinders have no
        // rotation)
        return {math::cos(phi), math::sin(phi), 0.f};
    }

};  // struct concentric_cylindrical2D

}  // namespace detray
