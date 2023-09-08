/** Detray library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/polar2D.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/propagator/detail/jacobian_kernel.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/tracks/detail/track_helper.hpp"
#include "detray/utils/invalid_values.hpp"

namespace detray::detail {

/// @brief Specialization for 2D cartesian frames
template <typename T, template <typename> class algebra_t>
struct jacobian_kernel<polar2D, T, algebra_t> {

    /// @name Type definitions for the struct
    /// @{
    using coordinate_frame = polar2D<algebra_t<T>>;
    using transform3_type = dtransform3D<algebra_t<T>>;
    using scalar_type = dscalar<algebra_t<T>>;
    using point2 = dpoint2D<algebra_t<T>>;
    using point3 = dpoint3D<algebra_t<T>>;
    using vector3 = dvector3D<algebra_t<T>>;
    // Matrix operator
    using matrix_operator = typename transform3_type::matrix_actor;
    // Matrix size type
    using size_type = typename matrix_operator::size_ty;
    // 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    // Rotation Matrix
    using rotation_matrix = matrix_type<3, 3>;
    // Shorthand vector/matrix types related to bound track parameters.
    using bound_vector = matrix_type<e_bound_size, 1>;
    using bound_matrix = matrix_type<e_bound_size, e_bound_size>;
    // Mapping from bound track parameters.
    using bound_to_free_matrix = matrix_type<e_free_size, e_bound_size>;
    // Shorthand vector/matrix types related to free track parameters.
    using free_vector = matrix_type<e_free_size, 1>;
    using free_matrix = matrix_type<e_free_size, e_free_size>;
    // Mapping from free track parameters.
    using free_to_bound_matrix = matrix_type<e_bound_size, e_free_size>;
    using free_to_path_matrix = matrix_type<1, e_free_size>;
    using path_to_free_matrix = matrix_type<e_free_size, 1>;
    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;
    /// @}

    DETRAY_HOST_DEVICE
    static inline rotation_matrix reference_frame(const transform3_type &trf3,
                                                  const point3 & /*pos*/,
                                                  const vector3 & /*dir*/) {
        return trf3.rotation();
    }

    DETRAY_HOST_DEVICE
    static inline free_to_path_matrix path_derivative(
        const transform3_type &trf3, const point3 & /*pos*/,
        const vector3 &dir) {

        free_to_path_matrix derivative =
            matrix_operator().template zero<1u, e_free_size>();

        const vector3 normal = coordinate_frame::normal(trf3);

        const vector3 pos_term = -1.f / vector::dot(normal, dir) * normal;

        matrix_operator().element(derivative, 0u, e_free_pos0) = pos_term[0];
        matrix_operator().element(derivative, 0u, e_free_pos1) = pos_term[1];
        matrix_operator().element(derivative, 0u, e_free_pos2) = pos_term[2];

        return derivative;
    }

    DETRAY_HOST_DEVICE
    static inline void set_bound_pos_to_free_pos_derivative(
        bound_to_free_matrix &bound_to_free_jacobian,
        const transform3_type &trf3, const point3 &pos, const vector3 &dir) {

        matrix_type<3, 2> bound_pos_to_free_pos_derivative =
            matrix_operator().template zero<3, 2>();

        const point3 local = coordinate_frame::global_to_local(trf3, pos, dir);
        const scalar_type lrad{local[0]};
        const scalar_type lphi{local[1]};

        const scalar_type lcos_phi{math_ns::cos(lphi)};
        const scalar_type lsin_phi{math_ns::sin(lphi)};

        // reference matrix
        const auto frame = reference_frame(trf3, pos, dir);

        // dxdu = d(x,y,z)/du
        const matrix_type<3, 1> dxdL =
            matrix_operator().template block<3, 1>(frame, 0u, 0u);
        // dxdv = d(x,y,z)/dv
        const matrix_type<3, 1> dydL =
            matrix_operator().template block<3, 1>(frame, 0u, 1u);

        const matrix_type<3, 1> col0 = dxdL * lcos_phi + dydL * lsin_phi;
        const matrix_type<3, 1> col1 =
            (dydL * lcos_phi - dxdL * lsin_phi) * lrad;

        matrix_operator().template set_block<3, 1>(
            bound_pos_to_free_pos_derivative, col0, e_free_pos0, e_bound_loc0);
        matrix_operator().template set_block<3, 1>(
            bound_pos_to_free_pos_derivative, col1, e_free_pos0, e_bound_loc1);

        matrix_operator().template set_block(bound_to_free_jacobian,
                                             bound_pos_to_free_pos_derivative,
                                             e_free_pos0, e_bound_loc0);
    }

    DETRAY_HOST_DEVICE
    static inline void set_free_pos_to_bound_pos_derivative(
        free_to_bound_matrix &free_to_bound_jacobian,
        const transform3_type &trf3, const point3 &pos, const vector3 &dir) {

        matrix_type<2, 3> free_pos_to_bound_pos_derivative =
            matrix_operator().template zero<2, 3>();

        const point3 local = coordinate_frame::global_to_local(trf3, pos, dir);

        const scalar_type lrad{local[0]};
        const scalar_type lphi{local[1]};

        const scalar_type lcos_phi{math_ns::cos(lphi)};
        const scalar_type lsin_phi{math_ns::sin(lphi)};

        // reference matrix
        const auto frame = reference_frame(trf3, pos, dir);
        const auto frameT = matrix_operator().transpose(frame);

        // dudG = du/d(x,y,z)
        const matrix_type<1, 3> dudG =
            matrix_operator().template block<1, 3>(frameT, 0u, 0u);
        // dvdG = dv/d(x,y,z)
        const matrix_type<1, 3> dvdG =
            matrix_operator().template block<1, 3>(frameT, 1u, 0u);

        const matrix_type<1, 3> row0 = dudG * lcos_phi + dvdG * lsin_phi;
        const matrix_type<1, 3> row1 =
            1. / lrad * (lcos_phi * dvdG - lsin_phi * dudG);

        matrix_operator().template set_block<1, 3>(
            free_pos_to_bound_pos_derivative, row0, e_bound_loc0, e_free_pos0);
        matrix_operator().template set_block<1, 3>(
            free_pos_to_bound_pos_derivative, row1, e_bound_loc1, e_free_pos0);

        matrix_operator().template set_block(free_to_bound_jacobian,
                                             free_pos_to_bound_pos_derivative,
                                             e_bound_loc0, e_free_pos0);
    }

    DETRAY_HOST_DEVICE
    static inline void set_bound_angle_to_free_pos_derivative(
        bound_to_free_matrix & /*bound_to_free_jacobian*/,
        const transform3_type & /*trf3*/, const point3 & /*pos*/,
        const vector3 & /*dir*/) {
        // Do nothing
    }
};

}  // namespace detray::detail
