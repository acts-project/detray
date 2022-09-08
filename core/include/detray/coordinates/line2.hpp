/** Detray library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinate_base.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

template <typename transform3_t>
struct line2 : public coordinate_base<line2, transform3_t> {

    /// @name Type definitions for the struct
    /// @{

    /// Base type
    using base_type = coordinate_base<line2, transform3_t>;
    /// Sclar type
    using scalar_type = typename base_type::scalar_type;
    /// Point in 2D space
    using point2 = typename base_type::point2;
    /// Point in 3D space
    using point3 = typename base_type::point3;
    /// Vector in 3D space
    using vector3 = typename base_type::vector3;

    /// @}

    /** This method transform from a point from 3D cartesian frame to a 2D
     * line point */
    DETRAY_HOST_DEVICE
    inline point2 operator()(const point3 &local3, scalar_type sign) const {

        return {sign * getter::perp(local3), local3[2]};
    }

    /** This method transform from a point from global cartesian 3D frame to a
     * local 2D line point */
    DETRAY_HOST_DEVICE
    inline point2 global_to_local(const transform3_t &trf, const point3 &p,
                                  const vector3 &d) const {

        const auto local3 = trf.point_to_local(p);

        // Line direction
        const vector3 z = trf.z();

        // Line center
        const point3 t = trf.translation();

        // Radial vector
        const vector3 r = vector::cross(z, d);

        // Assign the sign depending on the position w.r.t line
        // Right: -1
        // Left: 1
        const scalar_type sign =
            vector::dot(r, t - p) > 0. ? scalar_type{-1.} : scalar_type{1.};

        return this->operator()(local3, sign);
    }

    /** This method transform from a local 2D line point to a point global
     * cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 local_to_global(const transform3_t &trf,
                                                     const mask_t & /*mask*/,
                                                     const point2 &p,
                                                     const vector3 &d) const {

        // Line direction
        const vector3 z = trf.z();

        // Radial vector
        const vector3 r = vector::cross(z, d);

        // Local Z poisition in global cartesian coordinate
        const point3 locZ_in_global = trf.point_to_global(point3{0., 0., p[1]});

        return locZ_in_global + p[0] * vector::normalize(r);
    }
};

}  // namespace detray