/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"

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
    /// Point in 2D space
    using point2 = array_type<2>;

    /// Track vector types
    using bound_vector = matrix_type<e_bound_size, 1>;

    DETRAY_HOST_DEVICE
    inline void set_qop(bound_vector& bound_vec, const scalar_type& qop) {
        matrix_operator().element(bound_vec, e_bound_qoverp, 0u) = qop;
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
        const scalar_type sinTheta{math::sin(theta)};

        return {math::cos(phi) * sinTheta, math::sin(phi) * sinTheta,
                math::cos(theta)};
    }

    DETRAY_HOST_DEVICE
    inline scalar_type p(const bound_vector& bound_vec,
                         const scalar_type q) const {
        assert(qop(bound_vec) != 0.f);
        assert(q * qop(bound_vec) > 0.f);
        return q / qop(bound_vec);
    }

    DETRAY_HOST_DEVICE
    inline vector3 mom(const bound_vector& bound_vec,
                       const scalar_type q) const {
        return p(bound_vec, q) * dir(bound_vec);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qop(const bound_vector& bound_vec) const {
        return matrix_operator().element(bound_vec, e_bound_qoverp, 0u);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qopT(const bound_vector& bound_vec) const {
        const scalar_type theta{
            matrix_operator().element(bound_vec, e_bound_theta, 0u)};
        const scalar_type sinTheta{math::sin(theta)};
        assert(sinTheta != 0.f);
        return matrix_operator().element(bound_vec, e_bound_qoverp, 0u) /
               sinTheta;
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qopz(const bound_vector& bound_vec) const {
        const scalar_type theta{
            matrix_operator().element(bound_vec, e_bound_theta, 0u)};
        const scalar_type cosTheta{math::cos(theta)};
        assert(cosTheta != 0.f);
        return matrix_operator().element(bound_vec, e_bound_qoverp, 0u) /
               cosTheta;
    }

    DETRAY_HOST_DEVICE
    inline scalar_type time(const bound_vector& bound_vec) const {
        return matrix_operator().element(bound_vec, e_bound_time, 0u);
    }
};

}  // namespace detray::detail
