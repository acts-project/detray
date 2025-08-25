/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/geometry/coordinates/cylindrical2D.hpp"
#include "detray/propagator/detail/jacobian.hpp"

namespace detray::detail {

/// @brief Specialization for 2D cylindrical frames
template <concepts::algebra algebra_t>
struct jacobian<cylindrical2D<algebra_t>> {

    /// @name Type definitions for the struct
    /// @{
    using coordinate_frame = cylindrical2D<algebra_t>;

    using algebra_type = algebra_t;
    using transform3_type = dtransform3D<algebra_t>;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;

    // Rotation Matrix
    template <std::size_t ROWS, std::size_t COLS>
    using matrix_type = dmatrix<algebra_t, ROWS, COLS>;
    using rotation_matrix = dmatrix<algebra_t, 3, 3>;
    using bound_to_free_matrix_type = bound_to_free_matrix<algebra_t>;
    using free_to_bound_matrix_type = free_to_bound_matrix<algebra_t>;
    using free_to_path_matrix_type = free_to_path_matrix<algebra_t>;

    DETRAY_HOST_DEVICE
    static inline rotation_matrix reference_frame(const transform3_type &trf3,
                                                  const point3_type &pos,
                                                  const vector3_type &dir) {

        rotation_matrix rot = matrix::zero<rotation_matrix>();

        // y axis of the new frame is the z axis of cylindrical coordinate
        const vector3_type new_yaxis = trf3.z();

        // z axis of the new frame is the vector normal to the cylinder surface
        const point3_type local =
            coordinate_frame::global_to_local_3D(trf3, pos, dir);
        const vector3_type new_zaxis = coordinate_frame::normal(trf3, local);

        // x axis
        const vector3_type new_xaxis = vector::cross(new_yaxis, new_zaxis);

        getter::element(rot, 0u, 0u) = new_xaxis[0];
        getter::element(rot, 1u, 0u) = new_xaxis[1];
        getter::element(rot, 2u, 0u) = new_xaxis[2];
        getter::set_block(rot, new_yaxis, 0u, 1u);
        getter::element(rot, 0u, 2u) = new_zaxis[0];
        getter::element(rot, 1u, 2u) = new_zaxis[1];
        getter::element(rot, 2u, 2u) = new_zaxis[2];

        return rot;
    }

    DETRAY_HOST_DEVICE static inline free_to_path_matrix_type path_derivative(
        const transform3_type &trf3, const point3_type &pos,
        const vector3_type &dir, const vector3_type & /*dtds*/) {

        free_to_path_matrix_type derivative =
            matrix::zero<free_to_path_matrix_type>();

        const point3_type local =
            coordinate_frame::global_to_local_3D(trf3, pos, dir);
        const vector3_type normal = coordinate_frame::normal(trf3, local);

        const vector3_type pos_term =
            (-1.f / vector::dot(normal, dir)) * normal;

        getter::element(derivative, 0u, e_free_pos0) = pos_term[0];
        getter::element(derivative, 0u, e_free_pos1) = pos_term[1];
        getter::element(derivative, 0u, e_free_pos2) = pos_term[2];

        return derivative;
    }

    DETRAY_HOST_DEVICE
    static inline void set_bound_pos_to_free_pos_derivative(
        bound_to_free_matrix_type &bound_to_free_jacobian,
        const transform3_type &trf3, const point3_type &pos,
        const vector3_type &dir) {

        const rotation_matrix frame = reference_frame(trf3, pos, dir);

        // Get d(x,y,z)/d(loc0, loc1)
        const matrix_type<3, 2> bound_pos_to_free_pos_derivative =
            getter::block<3, 2>(frame, 0u, 0u);

        getter::set_block(bound_to_free_jacobian,
                          bound_pos_to_free_pos_derivative, e_free_pos0,
                          e_bound_loc0);
    }

    DETRAY_HOST_DEVICE
    static inline void set_free_pos_to_bound_pos_derivative(
        free_to_bound_matrix_type &free_to_bound_jacobian,
        const transform3_type &trf3, const point3_type &pos,
        const vector3_type &dir) {

        const rotation_matrix frame = reference_frame(trf3, pos, dir);
        const rotation_matrix frameT = matrix::transpose(frame);

        // Get d(loc0, loc1)/d(x,y,z)
        const matrix_type<2, 3> free_pos_to_bound_pos_derivative =
            getter::block<2, 3>(frameT, 0u, 0u);

        getter::set_block(free_to_bound_jacobian,
                          free_pos_to_bound_pos_derivative, e_bound_loc0,
                          e_free_pos0);
    }

    DETRAY_HOST_DEVICE
    static inline void set_bound_angle_to_free_pos_derivative(
        bound_to_free_matrix_type & /*bound_to_free_jacobian*/,
        const transform3_type & /*trf3*/, const point3_type & /*pos*/,
        const vector3_type & /*dir*/) {
        // Do nothing
    }
};

}  // namespace detray::detail
