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
#include "detray/definitions/track_parametrization.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

namespace detray {

/** Frame projection into a cartesian coordinate frame
 */
template <typename transform3_t>
struct cartesian2 final : public coordinate_base<cartesian2, transform3_t> {

    /// @name Type definitions for the struct
    /// @{

    // Transform type
    using transform3_type = transform3_t;
    // Base type
    using base_type = coordinate_base<cartesian2, transform3_t>;
    // Sclar type
    using scalar_type = typename base_type::scalar_type;
    // Point in 2D space
    using point2 = typename base_type::point2;
    // Point in 3D space
    using point3 = typename base_type::point3;
    // Vector in 3D space
    using vector3 = typename base_type::vector3;
    // Matrix operator
    using matrix_operator = typename base_type::matrix_operator;
    // Rotation Matrix
    using rotation_matrix = typename base_type::rotation_matrix;
    // Matrix size type
    using size_type = typename base_type::size_type;
    // 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type = typename base_type::template matrix_type<ROWS, COLS>;
    // Vector types
    using bound_vector = typename base_type::bound_vector;
    using free_vector = typename base_type::free_vector;
    // Matrix types
    using free_to_bound_matrix = typename base_type::free_to_bound_matrix;
    using bound_to_free_matrix = typename base_type::bound_to_free_matrix;
    using free_to_path_matrix = typename base_type::free_to_path_matrix;

    // Local point type in 2D cartesian coordinates
    using loc_point = point2;

    /// @}

    /** This method transforms a point from a global cartesian 3D frame to a
     * local 3D cartesian point */
    DETRAY_HOST_DEVICE
    inline point3 global_to_local(const transform3_t &trf3, const point3 &p,
                                  const vector3 & /*d*/) const {
        return trf3.point_to_local(p);
    }

    /** This method transforms a point from a global cartesian 3D frame to a
     * bound 2D cartesian point */
    DETRAY_HOST_DEVICE
    inline loc_point project_to_axes(const transform3_t &trf3, const point3 &p,
                                     const vector3 & /*d*/) const {
        auto loc_p = trf3.point_to_local(p);
        return {loc_p[0], loc_p[1]};
    }

    /** This method transform from a local 2D cartesian point to a point global
     * cartesian 3D frame*/
    DETRAY_HOST_DEVICE inline point3 local_to_global(const transform3_t &trf3,
                                                     const point3 &p) const {
        return trf3.point_to_global(p);
    }

    /** This method transform from a local 2D cartesian point to a point global
     * cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 bound_local_to_global(
        const transform3_t &trf3, const mask_t & /*mask*/, const point2 &p,
        const vector3 & /*d*/) const {

        return this->local_to_global(trf3, {p[0], p[1], 0.f});
    }

    /// @returns the normal vector
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline vector3 normal(const transform3_t &trf3,
                                             const point2 & = {},
                                             const mask_t & = {}) const {
        return trf3.z();
    }

    /// @returns the normal vector
    DETRAY_HOST_DEVICE inline vector3 normal(const transform3_t &trf3,
                                             const point3 & = {}) const {
        return trf3.z();
    }

    DETRAY_HOST_DEVICE inline rotation_matrix reference_frame(
        const transform3_t &trf3, const point3 & /*pos*/,
        const vector3 & /*dir*/) const {
        return trf3.rotation();
    }

    DETRAY_HOST_DEVICE inline free_to_path_matrix path_derivative(
        const transform3_t &trf3, const point3 & /*pos*/, const vector3 &dir,
        const vector3 & /*dtds*/) const {

        free_to_path_matrix derivative =
            matrix_operator().template zero<1u, e_free_size>();

        const vector3 normal = this->normal(trf3);

        const vector3 pos_term = -1.f / vector::dot(normal, dir) * normal;

        matrix_operator().element(derivative, 0u, e_free_pos0) = pos_term[0];
        matrix_operator().element(derivative, 0u, e_free_pos1) = pos_term[1];
        matrix_operator().element(derivative, 0u, e_free_pos2) = pos_term[2];

        return derivative;
    }

    DETRAY_HOST_DEVICE inline void set_bound_pos_to_free_pos_derivative(
        bound_to_free_matrix &bound_to_free_jacobian, const transform3_t &trf3,
        const point3 &pos, const vector3 &dir) const {

        const rotation_matrix frame = reference_frame(trf3, pos, dir);

        // Get d(x,y,z)/d(loc0, loc1)
        const matrix_type<3, 2> bound_pos_to_free_pos_derivative =
            matrix_operator().template block<3, 2>(frame, 0u, 0u);

        matrix_operator().template set_block(bound_to_free_jacobian,
                                             bound_pos_to_free_pos_derivative,
                                             e_free_pos0, e_bound_loc0);
    }

    DETRAY_HOST_DEVICE inline void set_free_pos_to_bound_pos_derivative(
        free_to_bound_matrix &free_to_bound_jacobian, const transform3_t &trf3,
        const point3 &pos, const vector3 &dir) const {

        const rotation_matrix frame = reference_frame(trf3, pos, dir);
        const rotation_matrix frameT = matrix_operator().transpose(frame);

        // Get d(loc0, loc1)/d(x,y,z)
        const matrix_type<2, 3> free_pos_to_bound_pos_derivative =
            matrix_operator().template block<2, 3>(frameT, 0, 0);

        matrix_operator().template set_block(free_to_bound_jacobian,
                                             free_pos_to_bound_pos_derivative,
                                             e_bound_loc0, e_free_pos0);
    }

    DETRAY_HOST_DEVICE inline void set_bound_angle_to_free_pos_derivative(
        bound_to_free_matrix & /*bound_to_free_jacobian*/,
        const transform3_t & /*trf3*/, const point3 & /*pos*/,
        const vector3 & /*dir*/) const {
        // Do nothing
    }
};  // struct cartesian2

}  // namespace detray
