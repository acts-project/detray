/** Detray plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {

/// @TODO: Move to algebra plugins
template <typename algebra_t>
struct matrix_helper {

    /// Matrix index type
    using size_type = dsize_type<algebra_t>;
    using scalar_type = dscalar<algebra_t>;
    using vector3 = dvector3D<algebra_t>;

    /// Matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type = dmatrix<algebra_t, ROWS, COLS>;

    /// Column-wise cross product between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> column_wise_cross(const matrix_type<3, 3>& m,
                                               const vector3& v) const {
        matrix_type<3, 3> ret;

        vector3 m_col0 = getter::vector<3>(m, 0, 0);
        vector3 m_col1 = getter::vector<3>(m, 0, 1);
        vector3 m_col2 = getter::vector<3>(m, 0, 2);

        getter::set_block(ret, static_cast<vector3>(vector::cross(m_col0, v)),
                          0, 0);
        getter::set_block(ret, static_cast<vector3>(vector::cross(m_col1, v)),
                          0, 1);
        getter::set_block(ret, static_cast<vector3>(vector::cross(m_col2, v)),
                          0, 2);

        return ret;
    }

    /// Column-wise multiplication between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> column_wise_multiply(const matrix_type<3, 3>& m,
                                                  const vector3& v) const {
        matrix_type<3, 3> ret;

        for (size_type i = 0; i < 3; i++) {
            for (size_type j = 0; j < 3; j++) {
                getter::element(ret, j, i) = getter::element(m, j, i) * v[j];
            }
        }

        return ret;
    }

    /// Cross product matrix
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> cross_matrix(const vector3& v) const {
        matrix_type<3, 3> ret;
        getter::element(ret, 0, 0) = 0;
        getter::element(ret, 0, 1) = -v[2];
        getter::element(ret, 0, 2) = v[1];
        getter::element(ret, 1, 0) = v[2];
        getter::element(ret, 1, 1) = 0;
        getter::element(ret, 1, 2) = -v[0];
        getter::element(ret, 2, 0) = -v[1];
        getter::element(ret, 2, 1) = v[0];
        getter::element(ret, 2, 2) = 0;

        return ret;
    }

    /// Outer product operation
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> outer_product(const vector3& v1,
                                           const vector3& v2) const {
        matrix_type<3, 1> m1;
        getter::element(m1, 0, 0) = v1[0];
        getter::element(m1, 1, 0) = v1[1];
        getter::element(m1, 2, 0) = v1[2];

        matrix_type<1, 3> m2;
        getter::element(m2, 0, 0) = v2[0];
        getter::element(m2, 0, 1) = v2[1];
        getter::element(m2, 0, 2) = v2[2];

        return m1 * m2;
    }

    /// Cholesky decompose
    template <size_type N>
    DETRAY_HOST_DEVICE inline matrix_type<N, N> cholesky_decompose(
        const matrix_type<N, N>& mat) const {

        matrix_type<N, N> L = matrix::zero<matrix_type<N, N>>();

        // Cholesky–Banachiewicz algorithm
        for (size_type i = 0u; i < N; i++) {
            for (size_type j = 0u; j <= i; j++) {
                scalar_type sum = 0.f;
                for (size_type k = 0u; k < j; k++)
                    sum += getter::element(L, i, k) * getter::element(L, j, k);

                if (i == j) {
                    getter::element(L, i, j) = static_cast<scalar_type>(
                        math::sqrt(getter::element(mat, i, i) - sum));
                } else {
                    getter::element(L, i, j) =
                        (1.f / getter::element(L, j, j) *
                         (getter::element(mat, i, j) - sum));
                }
            }
        }

        return L;
    }
};

}  // namespace detray
