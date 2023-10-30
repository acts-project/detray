/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"

// System include(s).
#include <cmath>

namespace detray::detail {

template <typename matrix_operator_t>
struct track_helper {

    /// Matrix actor
    using matrix_operator = matrix_operator_t;
    /// Size type
    using size_type = typename matrix_operator_t::size_ty;
    /// Scalar type
    using scalar_type = typename matrix_operator_t::scalar_type;
    /// 2D Matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    /// Array type
    template <size_type N>
    using array_type = typename matrix_operator::template array_type<N>;
    /// 3-element "vector" type
    using vector3 = array_type<3>;
    /// Point in 3D space
    using point3 = vector3;
    /// Point in 2D space
    using point2 = array_type<2>;

    /// Track vector types
    using bound_vector = matrix_type<e_bound_size, 1>;
    using free_vector = matrix_type<e_free_size, 1>;

    DETRAY_HOST_DEVICE
    inline point3 pos(const free_vector& free_vec) const {
        return {matrix_operator().element(free_vec, e_free_pos0, 0u),
                matrix_operator().element(free_vec, e_free_pos1, 0u),
                matrix_operator().element(free_vec, e_free_pos2, 0u)};
    }

    DETRAY_HOST_DEVICE
    inline void set_pos(free_vector& free_vec, const point3& pos) {
        matrix_operator().element(free_vec, e_free_pos0, 0u) = pos[0];
        matrix_operator().element(free_vec, e_free_pos1, 0u) = pos[1];
        matrix_operator().element(free_vec, e_free_pos2, 0u) = pos[2];
    }

    DETRAY_HOST_DEVICE
    inline vector3 dir(const free_vector& free_vec) const {
        return {matrix_operator().element(free_vec, e_free_dir0, 0u),
                matrix_operator().element(free_vec, e_free_dir1, 0u),
                matrix_operator().element(free_vec, e_free_dir2, 0u)};
    }

    DETRAY_HOST_DEVICE
    inline void set_dir(free_vector& free_vec, const vector3& dir) {
        matrix_operator().element(free_vec, e_free_dir0, 0u) = dir[0];
        matrix_operator().element(free_vec, e_free_dir1, 0u) = dir[1];
        matrix_operator().element(free_vec, e_free_dir2, 0u) = dir[2];
    }

    DETRAY_HOST_DEVICE
    inline point2 bound_local(const bound_vector& bound_vec) const {
        return {matrix_operator().element(bound_vec, e_bound_loc0, 0u),
                matrix_operator().element(bound_vec, e_bound_loc1, 0u)};
    }

    DETRAY_HOST_DEVICE
    inline vector3 dir(const bound_vector& bound_vec) const {
        const scalar_type phi{
            matrix_operator().element(bound_vec, e_bound_phi, 0u)};
        const scalar_type theta{
            matrix_operator().element(bound_vec, e_bound_theta, 0u)};
        const scalar_type sinTheta{math_ns::sin(theta)};

        return {math_ns::cos(phi) * sinTheta, math_ns::sin(phi) * sinTheta,
                math_ns::cos(theta)};
    }

    DETRAY_HOST_DEVICE
    inline scalar_type p(const free_vector& free_vec) const {
        return charge(free_vec) / qop(free_vec);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type p(const bound_vector& bound_vec) const {
        return charge(bound_vec) / qop(bound_vec);
    }

    DETRAY_HOST_DEVICE
    inline vector3 mom(const free_vector& free_vec) const {
        return p(free_vec) * dir(free_vec);
    }

    DETRAY_HOST_DEVICE
    inline vector3 mom(const bound_vector& bound_vec) const {
        return p(bound_vec) * dir(bound_vec);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qop(const free_vector& free_vec) const {
        return matrix_operator().element(free_vec, e_free_qoverp, 0u);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qop(const bound_vector& bound_vec) const {
        return matrix_operator().element(bound_vec, e_bound_qoverp, 0u);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qopT(const free_vector& free_vec) const {
        const auto dir = this->dir(free_vec);
        return matrix_operator().element(free_vec, e_free_qoverp, 0u) /
               getter::perp(dir);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qopT(const bound_vector& bound_vec) const {
        const scalar_type theta{
            matrix_operator().element(bound_vec, e_bound_theta, 0u)};
        const scalar_type sinTheta{math_ns::sin(theta)};
        return matrix_operator().element(bound_vec, e_bound_qoverp, 0u) /
               sinTheta;
    }

    DETRAY_HOST_DEVICE
    inline scalar_type time(const free_vector& free_vec) const {
        return matrix_operator().element(free_vec, e_free_time, 0u);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type time(const bound_vector& bound_vec) const {
        return matrix_operator().element(bound_vec, e_bound_time, 0u);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type charge(const free_vector& free_vec) const {
        return std::signbit(
                   matrix_operator().element(free_vec, e_free_qoverp, 0u))
                   ? -1.f
                   : 1.f;
    }

    DETRAY_HOST_DEVICE
    inline scalar_type charge(const bound_vector& bound_vec) const {
        return std::signbit(
                   matrix_operator().element(bound_vec, e_bound_qoverp, 0u))
                   ? -1.f
                   : 1.f;
    }
};

}  // namespace detray::detail
