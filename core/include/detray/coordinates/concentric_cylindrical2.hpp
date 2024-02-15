/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinate_base.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/utils/invalid_values.hpp"

namespace detray {

/// Local frame projection into a 2D concentric cylindrical coordinate frame
/// (No rotation in coordinate transformation)
template <typename transform3_t>
struct concentric_cylindrical2
    : public coordinate_base<concentric_cylindrical2, transform3_t> {

    /// @name Type definitions for the struct
    /// @{

    // Transform type
    using transform3_type = transform3_t;
    // Base type
    using base_type = coordinate_base<concentric_cylindrical2, transform3_t>;
    // Sclar type
    using scalar_type = typename transform3_t::scalar_type;
    // Point in 2D space
    using point2 = typename transform3_t::point2;
    // Point in 3D space
    using point3 = typename transform3_t::point3;
    // Vector in 3D space
    using vector3 = typename transform3_t::vector3;
    // Matrix actor
    using matrix_operator = typename base_type::matrix_operator;
    // Matrix size type
    using size_type = typename base_type::size_type;
    // Rotation Matrix
    using rotation_matrix = typename base_type::rotation_matrix;
    // Vector types
    using bound_vector = typename base_type::bound_vector;
    using free_vector = typename base_type::free_vector;
    // Matrix types
    using free_to_bound_matrix = typename base_type::free_to_bound_matrix;
    using bound_to_free_matrix = typename base_type::bound_to_free_matrix;
    using free_to_path_matrix = typename base_type::free_to_path_matrix;

    // Local point type in 2D cylindrical coordinates
    using loc_point = point2;

    /// @}

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 2D cylindrical point
    DETRAY_HOST_DEVICE
    inline point3 global_to_local(const transform3_t &trf, const point3 &p,
                                  const vector3 & /*d*/) const {
        const point3 local3 = p - trf.translation();

        return {getter::perp(local3) * getter::phi(local3), local3[2],
                getter::perp(local3)};
    }

    /// This method transforms a point from a global cartesian 3D frame to a
    /// bound 2D cylindrical point
    DETRAY_HOST_DEVICE
    inline loc_point project_to_axes(const transform3_t &trf, const point3 &p,
                                     const vector3 & /*d*/) const {
        const point3 local3 = trf.point_to_local(p);

        return {getter::phi(local3), local3[2]};
    }

    /// This method transform from a local 2D cylindrical point to a point
    /// global cartesian 3D frame
    DETRAY_HOST_DEVICE inline point3 local_to_global(const transform3_t &trf,
                                                     const point3 &p) const {
        const scalar_type r{p[2]};
        const scalar_type phi{p[0] / r};
        const scalar_type x{r * math::cos(phi)};
        const scalar_type y{r * math::sin(phi)};
        const scalar_type z{p[1]};

        return point3{x, y, z} + trf.translation();
    }

    /// This method transform from a local 2D cylindrical point to a point
    /// global cartesian 3D frame
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 bound_local_to_global(
        const transform3_t &trf, const mask_t &mask, const point2 &p,
        const vector3 & /*dir*/) const {

        return this->local_to_global(trf,
                                     {p[0], p[1], mask[mask_t::shape::e_r]});
    }

    /// @returns the normal vector
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline vector3 normal(const transform3_t &,
                                             const point2 &bound_pos,
                                             const mask_t &mask) const {
        const scalar_type phi{bound_pos[0] / mask[mask_t::shape::e_r]};
        const vector3 local_normal{math::cos(phi), math::sin(phi), 0.f};

        // normal vector in local coordinate
        return local_normal;
    }

    /// @returns the normal vector given a local position @param loc_pos
    DETRAY_HOST_DEVICE inline vector3 normal(const transform3_t &,
                                             const point3 &loc_pos) const {
        const scalar_type phi{loc_pos[0] / loc_pos[2]};
        const vector3 local_normal{math::cos(phi), math::sin(phi), 0.f};

        // normal vector in local coordinate
        return local_normal;
    }

    DETRAY_HOST_DEVICE inline rotation_matrix reference_frame(
        const transform3_t &trf3, const point3 &pos, const vector3 &dir) const {

        rotation_matrix rot = matrix_operator().template zero<3, 3>();

        // y axis of the new frame is the z axis of cylindrical coordinate
        const auto new_yaxis = vector3{0.f, 0.f, 1.f};

        // z axis of the new frame is the vector normal to the cylinder surface
        const point3 local = this->global_to_local(trf3, pos, dir);
        const vector3 new_zaxis = normal(trf3, local);

        // x axis
        const vector3 new_xaxis = vector::cross(new_yaxis, new_zaxis);

        matrix_operator().element(rot, 0u, 0u) = new_xaxis[0];
        matrix_operator().element(rot, 1u, 0u) = new_xaxis[1];
        matrix_operator().element(rot, 2u, 0u) = new_xaxis[2];
        matrix_operator().element(rot, 0u, 1u) = new_yaxis[0];
        matrix_operator().element(rot, 1u, 1u) = new_yaxis[1];
        matrix_operator().element(rot, 2u, 1u) = new_yaxis[2];
        matrix_operator().element(rot, 0u, 2u) = new_zaxis[0];
        matrix_operator().element(rot, 1u, 2u) = new_zaxis[1];
        matrix_operator().element(rot, 2u, 2u) = new_zaxis[2];

        return rot;
    }

    DETRAY_HOST_DEVICE inline free_to_path_matrix path_derivative(
        const transform3_t &trf3, const point3 &pos, const vector3 &dir,
        const vector3 & /*dtds*/) const {

        free_to_path_matrix derivative =
            matrix_operator().template zero<1u, e_free_size>();

        const point3 local = this->global_to_local(trf3, pos, dir);
        const vector3 normal = this->normal(trf3, local);

        const vector3 pos_term = -1.f / vector::dot(normal, dir) * normal;

        matrix_operator().element(derivative, 0u, e_free_pos0) = pos_term[0];
        matrix_operator().element(derivative, 0u, e_free_pos1) = pos_term[1];
        matrix_operator().element(derivative, 0u, e_free_pos2) = pos_term[2];

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
        bound_to_free_matrix & /*bound_to_free_jacobian*/,
        const transform3_t & /*trf3*/, const point3 & /*pos*/,
        const vector3 & /*dir*/) const {
        // Do nothing
    }
};  // struct concentric_cylindrical2

}  // namespace detray
