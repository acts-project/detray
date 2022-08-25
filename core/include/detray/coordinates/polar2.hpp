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

template <typename transform3_t>
struct polar2 : public coordinate_base<polar2, transform3_t> {

    /// @name Type definitions for the struct
    /// @{

    // Base type
    using base_type = coordinate_base<polar2, transform3_t>;
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

    /** This method transform from a point from 2D cartesian frame to a 2D
     * polar point */
    DETRAY_HOST_DEVICE inline point2 operator()(const point2 &p) const {

        return {getter::perp(p), getter::phi(p)};
    }

    /** This method transform from a point in 3D cartesian frame to a 2D
     * polar point */
    DETRAY_HOST_DEVICE inline point2 operator()(const point3 &p) const {

        return {getter::perp(p), getter::phi(p)};
    }

    /** This method transform from a point from global cartesian 3D frame to a
     * local 2D polar point */
    DETRAY_HOST_DEVICE
    inline point2 global_to_local(const transform3_t &trf, const point3 &p,
                                  const vector3 & /*d*/) const {
        const auto local3 = trf.point_to_local(p);
        return this->operator()(local3);
    }

    /** This method transform from a local 2D polar point to a point global
     * cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 local_to_global(
        const transform3_t &trf, const mask_t & /*mask*/, const point2 &p,
        const vector3 & /*d*/) const {
        const scalar_type x = p[0] * std::cos(p[1]);
        const scalar_type y = p[0] * std::sin(p[1]);

        return trf.point_to_global(point3{x, y, 0.});
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline rotation_matrix reference_frame(
        const transform3_t &trf3, const mask_t & /*mask*/,
        const point3 & /*pos*/, const vector3 & /*dir*/) const {
        return trf3.rotation();
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline matrix_type<3, 2> bound_to_free_rotation(
        const transform3_t &trf3, const mask_t &mask, const point3 &pos,
        const vector3 &dir) const {

        matrix_type<3, 2> bound_to_free_rotation =
            matrix_actor().template zero<3, 2>();

        const point2 local2 = this->operator()(pos);
        const scalar_type lrad = local2[0];
        const scalar_type lphi = local2[1];

        const scalar_type lcos_phi = std::cos(lphi);
        const scalar_type lsin_phi = std::sin(lphi);

        // reference matrix
        const auto frame = reference_frame(trf3, pos, dir);

        // dxdL = dx/d(u,v,w)
        const auto dxdL = matrix_actor().template block<3, 1>(frame, 0, 0);
        // dydL = dy/d(u,v,w)
        const auto dydL = matrix_actor().template block<3, 1>(frame, 0, 1);

        const auto col0 = dxdL * lcos_phi + dydL * lsin_phi;
        const auto col1 = (dydL * lcos_phi - dxdL * lsin_phi) * lrad;

        matrix_actor().set_block<3, 1>(bound_to_free_rotation, col0,
                                       e_free_pos0, e_bound_loc0);
        matrix_actor().set_block<3, 1>(bound_to_free_rotation, col1,
                                       e_free_pos0, e_bound_loc1);

        return bound_to_free_rotation;
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline matrix_type<2, 3> free_to_bound_rotation(
        const transform3_t &trf3, const mask_t &mask, const point3 &pos,
        const vector3 &dir) const {

        matrix_type<2, 3> free_to_bound_rotation =
            matrix_actor().template zero<2, 3>();

        const auto local = this->global_to_local(trf3, pos, dir);

        const scalar_type lrad = local[0];
        const scalar_type lphi = local[1];

        const scalar_type lcos_phi = std::cos(lphi);
        const scalar_type lsin_phi = std::sin(lphi);

        // reference matrix
        const auto frame = reference_frame(trf3, pos, dir);
        const auto frameT = matrix_actor().transpose(frame);

        // dudG = du/d(x,y,z)
        const auto dudG = matrix_actor().template block<3, 1>(frameT, 0, 0);
        // dvdG = dv/d(x,y,z)
        const auto dvdG = matrix_actor().template block<3, 1>(frameT, 0, 1);

        const auto row0 = dudG * lcos_phi + dvdG * lsin_phi;
        const auto row1 = (dvdG * lcos_phi - dudG * lsin_phi) * 1. / lrad;

        matrix_actor().set_block<1, 3>(free_to_bound_rotation, row0,
                                       e_bound_loc0, e_free_pos0);
        matrix_actor().set_block<1, 3>(free_to_bound_rotation, row1,
                                       e_bound_loc1, e_free_pos0);

        return free_to_bound_rotation;
    }
};

}  // namespace detray