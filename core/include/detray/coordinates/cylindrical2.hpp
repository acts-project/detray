/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinate_base.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s).
#include <cmath>

namespace detray {

/** Local frame projection into a polar coordinate frame
 */
template <typename transform3_t>
struct cylindrical2 : public coordinate_base<cylindrical2, transform3_t> {

    /// @name Type definitions for the struct
    /// @{

    // Base type
    using base_type = coordinate_base<cylindrical2, transform3_t>;
    // Sclar type
    using scalar_type = typename transform3_t::scalar_type;
    // Point in 2D space
    using point2 = typename transform3_t::point2;
    // Point in 3D space
    using point3 = typename transform3_t::point3;
    // Vector in 3D space
    using vector3 = typename transform3_t::vector3;
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
    // Matrix types
    using free_to_bound_matrix = typename base_type::free_to_bound_matrix;
    using bound_to_free_matrix = typename base_type::bound_to_free_matrix;

    /// @}

    /** This method transform from a point from 3D cartesian frame to a 2D
     * cylindrical point */
    DETRAY_HOST_DEVICE
    inline point2 operator()(const point3 &p) const {

        return {getter::perp(p) * getter::phi(p), p[2]};
    }

    /** This method transform from a point from global cartesian 3D frame to a
     * local 2D cylindrical point */
    DETRAY_HOST_DEVICE
    inline point2 global_to_local(const transform3_t &trf, const point3 &p,
                                  const vector3 & /*d*/) const {
        const auto local3 = trf.point_to_local(p);
        return this->operator()(local3);
    }

    /** This method transform from a local 2D cylindrical point to a point
     * global cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 local_to_global(
        const transform3_t &trf, const mask_t &mask, const point2 &p,
        const vector3 & /*d*/) const {
        const scalar_type r = mask.radius();
        const scalar_type phi = p[0] / r;
        const scalar_type x = r * std::cos(phi);
        const scalar_type y = r * std::sin(phi);
        const scalar_type z = p[1];

        return trf.point_to_global(point3{x, y, z});
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline vector3 normal(const transform3_t &trf3,
                                             const mask_t &mask,
                                             const point3 &pos,
                                             const vector3 &dir) const {
        const point2 local2 = this->global_to_local(trf3, pos, dir);
        const scalar_type r = mask.radius();
        const scalar_type phi = local2[0] / r;

        // normal vector in local coordinate
        const vector3 local_normal{std::cos(phi), std::sin(phi), 0};

        // normal vector in global coordinate
        return trf3.rotation() * local_normal;
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline rotation_matrix reference_frame(
        const transform3_t &trf3, const mask_t &mask, const point3 &pos,
        const vector3 &dir) const {

        rotation_matrix rot = matrix_actor().template zero<3, 3>();

        // y axis of the new frame is the z axis of cylindrical coordinate
        const auto new_yaxis =
            matrix_actor().template block<3, 1>(trf3.matrix(), 0, 2);

        // z axis of the new frame is the vector normal to the cylinder surface
        const vector3 new_zaxis = normal(trf3, mask, pos, dir);

        // x axis
        const vector3 new_xaxis = vector::cross(new_yaxis, new_zaxis);

        matrix_actor().element(rot, 0, 0) = new_xaxis[0];
        matrix_actor().element(rot, 1, 0) = new_xaxis[1];
        matrix_actor().element(rot, 2, 0) = new_xaxis[2];
        matrix_actor().template set_block<3, 1>(rot, new_yaxis, 0, 1);
        matrix_actor().element(rot, 0, 2) = new_zaxis[0];
        matrix_actor().element(rot, 1, 2) = new_zaxis[1];
        matrix_actor().element(rot, 2, 2) = new_zaxis[2];

        return rot;
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline void set_bound_pos_to_free_pos_derivative(
        bound_to_free_matrix &free_to_bound_jacobian, const transform3_t &trf3,
        const mask_t &mask, const point3 &pos, const vector3 &dir) const {

        const auto frame = reference_frame(trf3, mask, pos, dir);

        // Get d(x,y,z)/d(loc0, loc1)
        const auto bound_pos_to_free_pos_derivative =
            matrix_actor().template block<3, 2>(frame, 0, 0);

        matrix_actor().template set_block(free_to_bound_jacobian,
                                          bound_pos_to_free_pos_derivative,
                                          e_free_pos0, e_bound_loc0);
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline void set_free_pos_to_bound_pos_derivative(
        free_to_bound_matrix &bound_to_free_jacobian, const transform3_t &trf3,
        const mask_t &mask, const point3 &pos, const vector3 &dir) const {

        const auto frame = reference_frame(trf3, mask, pos, dir);
        const auto frameT = matrix_actor().transpose(frame);

        // Get d(loc0, loc1)/d(x,y,z)
        const auto free_pos_to_bound_pos_derivative =
            matrix_actor().template block<2, 3>(frameT, 0, 0);

        matrix_actor().template set_block(bound_to_free_jacobian,
                                          free_pos_to_bound_pos_derivative,
                                          e_bound_loc0, e_free_pos0);
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline void set_bound_angle_to_free_pos_derivative(
        bound_to_free_matrix & /*bound_to_free_jacobian*/,
        const transform3_t & /*trf3*/, const mask_t & /*mask*/,
        const point3 & /*pos*/, const vector3 & /*dir*/) const {
        // Do nothing
    }
};  // struct cylindrical2

}  // namespace detray