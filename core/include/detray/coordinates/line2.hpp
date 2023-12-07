/** Detray library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinate_base.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

namespace detray {

template <typename transform3_t>
struct line2 : public coordinate_base<line2, transform3_t> {

    /// @name Type definitions for the struct
    /// @{

    // Transform type
    using transform3_type = transform3_t;
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
    using matrix_operator = typename base_type::matrix_operator;
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
    // Matrix types
    using free_to_bound_matrix = typename base_type::free_to_bound_matrix;
    using bound_to_free_matrix = typename base_type::bound_to_free_matrix;
    using free_to_path_matrix = typename base_type::free_to_path_matrix;

    // Local point type in line coordinates
    using loc_point = point2;

    /// @}

    /** This method transform from a point from global cartesian 3D frame to a
     * local 2D line point */
    DETRAY_HOST_DEVICE
    inline point3 global_to_local(const transform3_t &trf, const point3 &p,
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
        const scalar_type sign = vector::dot(r, t - p) > 0.f ? -1.f : 1.f;

        return {sign * getter::perp(local3), local3[2], getter::phi(local3)};
    }

    /** This method transform from a local 2D line point to a point global
     * cartesian 3D frame*/
    DETRAY_HOST_DEVICE inline point3 local_to_global(const transform3_t &trf,
                                                     const point3 &p) const {
        const scalar_type R = std::abs(p[0]);
        const point3 local = {R * math_ns::cos(p[2]), R * math_ns::sin(p[2]),
                              p[1]};

        return trf.point_to_global(local);
    }

    /** This method transform from a local 2D line point to a point global
     * cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 bound_local_to_global(
        const transform3_t &trf, const mask_t & /*mask*/, const point2 &p,
        const vector3 &d) const {

        // Line direction
        const vector3 z = trf.z();

        // Radial vector
        const vector3 r = vector::cross(z, d);

        // Local Z poisition in global cartesian coordinate
        const point3 locZ_in_global =
            trf.point_to_global(point3{0.f, 0.f, p[1]});

        return locZ_in_global + p[0] * vector::normalize(r);
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
        const vector3 &dir) const {

        rotation_matrix rot = matrix_operator().template zero<3, 3>();

        // y axis of the new frame is the z axis of line coordinate
        const auto new_yaxis =
            matrix_operator().template block<3, 1>(trf3.matrix(), 0u, 2u);

        // x axis of the new frame is (yaxis x track direction)
        const auto new_xaxis = vector::normalize(vector::cross(new_yaxis, dir));

        // z axis
        const auto new_zaxis = vector::cross(new_xaxis, new_yaxis);

        matrix_operator().element(rot, 0u, 0u) = new_xaxis[0];
        matrix_operator().element(rot, 1u, 0u) = new_xaxis[1];
        matrix_operator().element(rot, 2u, 0u) = new_xaxis[2];
        matrix_operator().template set_block<3, 1>(rot, new_yaxis, 0u, 1u);
        matrix_operator().element(rot, 0u, 2u) = new_zaxis[0];
        matrix_operator().element(rot, 1u, 2u) = new_zaxis[1];
        matrix_operator().element(rot, 2u, 2u) = new_zaxis[2];

        return rot;
    }

    DETRAY_HOST_DEVICE inline free_to_path_matrix path_derivative(
        const transform3_t &trf3, const point3 &pos, const vector3 &dir,
        const vector3 &dtds) const {

        free_to_path_matrix derivative =
            matrix_operator().template zero<1u, e_free_size>();

        // The vector between position and center
        const point3 center = trf3.translation();
        const vector3 pc = pos - center;

        // The local frame z axis
        const vector3 local_zaxis = getter::vector<3>(trf3.matrix(), 0u, 2u);

        // The local z coordinate
        const scalar_type pz = vector::dot(pc, local_zaxis);

        // Cosine of angle between momentum direction and local frame z axis
        const scalar_type dz = vector::dot(local_zaxis, dir);

        // local x axis component of pc:
        const vector3 pc_x = pc - pz * local_zaxis;

        const scalar_type norm =
            -1.f / (1.f - dz * dz + vector::dot(pc_x, dtds));

        const vector3 pos_term = norm * (dir - dz * local_zaxis);
        const vector3 dir_term = norm * pc_x;

        matrix_operator().element(derivative, 0u, e_free_pos0) = pos_term[0];
        matrix_operator().element(derivative, 0u, e_free_pos1) = pos_term[1];
        matrix_operator().element(derivative, 0u, e_free_pos2) = pos_term[2];
        matrix_operator().element(derivative, 0u, e_free_dir0) = dir_term[0];
        matrix_operator().element(derivative, 0u, e_free_dir1) = dir_term[1];
        matrix_operator().element(derivative, 0u, e_free_dir2) = dir_term[2];

        return derivative;
    }

    DETRAY_HOST_DEVICE inline void set_bound_pos_to_free_pos_derivative(
        bound_to_free_matrix &bound_to_free_jacobian, const transform3_t &trf3,
        const point3 &pos, const vector3 &dir) const {

        const auto frame = reference_frame(trf3, pos, dir);

        // Get d(x,y,z)/d(loc0, loc1)
        const auto bound_pos_to_free_pos_derivative =
            matrix_operator().template block<3, 2>(frame, 0u, 0u);

        matrix_operator().template set_block(bound_to_free_jacobian,
                                             bound_pos_to_free_pos_derivative,
                                             e_free_pos0, e_bound_loc0);
    }

    DETRAY_HOST_DEVICE inline void set_free_pos_to_bound_pos_derivative(
        free_to_bound_matrix &free_to_bound_jacobian, const transform3_t &trf3,
        const point3 &pos, const vector3 &dir) const {

        const auto frame = reference_frame(trf3, pos, dir);
        const auto frameT = matrix_operator().transpose(frame);

        // Get d(loc0, loc1)/d(x,y,z)
        const auto free_pos_to_bound_pos_derivative =
            matrix_operator().template block<2, 3>(frameT, 0u, 0u);

        matrix_operator().template set_block(free_to_bound_jacobian,
                                             free_pos_to_bound_pos_derivative,
                                             e_bound_loc0, e_free_pos0);
    }

    DETRAY_HOST_DEVICE inline void set_bound_angle_to_free_pos_derivative(
        bound_to_free_matrix &bound_to_free_jacobian, const transform3_t &trf3,
        const point3 &pos, const vector3 &dir) const {

        // local2
        const auto local2 = this->global_to_local(trf3, pos, dir);

        // Reference frame
        const auto frame = reference_frame(trf3, pos, dir);

        // new x_axis
        const auto new_xaxis = getter::vector<3>(frame, 0u, 0u);

        // new y_axis
        const auto new_yaxis = getter::vector<3>(frame, 0u, 1u);

        // new z_axis
        const auto new_zaxis = getter::vector<3>(frame, 0u, 2u);

        // the projection of direction onto ref frame normal
        const scalar_type ipdn{1.f / vector::dot(dir, new_zaxis)};

        // d(n_x,n_y,n_z)/dPhi
        const auto dNdPhi = matrix_operator().template block<3, 1>(
            bound_to_free_jacobian, e_free_dir0, e_bound_phi);

        // Get new_yaxis X d(n_x,n_y,n_z)/dPhi
        auto y_cross_dNdPhi = vector::cross(new_yaxis, dNdPhi);

        // d(n_x,n_y,n_z)/dTheta
        const auto dNdTheta = matrix_operator().template block<3, 1>(
            bound_to_free_jacobian, e_free_dir0, e_bound_theta);

        // build the cross product of d(D)/d(eBoundPhi) components with y axis
        auto y_cross_dNdTheta = vector::cross(new_yaxis, dNdTheta);

        const scalar_type C{ipdn * local2[0]};
        // and correct for the x axis components
        vector3 phi_to_free_pos_derivative =
            y_cross_dNdPhi - new_xaxis * vector::dot(new_xaxis, y_cross_dNdPhi);

        phi_to_free_pos_derivative = C * phi_to_free_pos_derivative;

        vector3 theta_to_free_pos_derivative =
            y_cross_dNdTheta -
            new_xaxis * vector::dot(new_xaxis, y_cross_dNdTheta);

        theta_to_free_pos_derivative = C * theta_to_free_pos_derivative;

        // Set the jacobian components
        matrix_operator().element(bound_to_free_jacobian, e_free_pos0,
                                  e_bound_phi) = phi_to_free_pos_derivative[0];
        matrix_operator().element(bound_to_free_jacobian, e_free_pos1,
                                  e_bound_phi) = phi_to_free_pos_derivative[1];
        matrix_operator().element(bound_to_free_jacobian, e_free_pos2,
                                  e_bound_phi) = phi_to_free_pos_derivative[2];
        matrix_operator().element(bound_to_free_jacobian, e_free_pos0,
                                  e_bound_theta) =
            theta_to_free_pos_derivative[0];
        matrix_operator().element(bound_to_free_jacobian, e_free_pos1,
                                  e_bound_theta) =
            theta_to_free_pos_derivative[1];
        matrix_operator().element(bound_to_free_jacobian, e_free_pos2,
                                  e_bound_theta) =
            theta_to_free_pos_derivative[2];
    }
};

}  // namespace detray
