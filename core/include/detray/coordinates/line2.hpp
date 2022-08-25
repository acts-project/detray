/** Algebra plugins library, part of the ACTS project
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

    // Base type
    using base_type = coordinate_base<line2, transform3_t>;
    // Sclar type
    using scalar_type = typename base_type::scalar_type;
    // Point in 2D space
    using point2 = typename base_type::point2;
    // Point in 3D space
    using point3 = typename base_type::point3;
    // Vector in 3D space
    using vector3 = typename base_type::vector3;
    // Matrix actor
    using matrix_actor = typename base_type::matrix_actor;
    // Matrix size type
    using size_type = typename base_type::size_type;
    // 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type = typename base_type::template matrix_type<ROWS, COLS>;
    // Rotation Matrix
    using rotation_matrix = typename base_type::rotation_matrix;
    // Vector types
    using bound_vector = typename base_type::bound_vector;
    using free_vector = typename base_type::free_vector;
    // Track Helper
    using track_helper = typename base_type::track_helper;

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

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline rotation_matrix reference_frame(
        const transform3_t &trf3, const mask_t & /*mask*/,
        const point3 & /*pos*/, const vector3 &dir) const {

        rotation_matrix rot = matrix_actor().template zero<3, 3>();

        // y axis of the new frame is the z axis of line coordinate
        const auto new_yaxis = matrix_actor().template block<3, 1>(trf3, 0, 2);

        // x axis of the new frame is (yaxis x track direction)
        const auto new_xaxis = vector::cross(new_yaxis, dir);

        // z axis
        const auto new_zaxis = vector::cross(new_xaxis, new_yaxis);

        matrix_actor().set_block<3, 1>(rot, new_xaxis, 0, 0);
        matrix_actor().set_block<3, 1>(rot, new_yaxis, 0, 1);
        matrix_actor().set_block<3, 1>(rot, new_zaxis, 0, 2);

        return rot;
    }
};

}  // namespace detray