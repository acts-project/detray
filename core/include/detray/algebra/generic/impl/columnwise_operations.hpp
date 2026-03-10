/** Detray plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"

namespace detray {

/// @TODO: Move to algebra plugins
template <concepts::algebra algebra_t>
struct matrix_helper {

    /// Matrix index type
    using index_type = dindex_type<algebra_t>;
    using scalar_type = dscalar<algebra_t>;
    using vector3 = dvector3D<algebra_t>;

    /// Matrix type
    template <index_type ROWS, index_type COLS>
    using matrix_type = dmatrix<algebra_t, ROWS, COLS>;

    /// Column-wise cross product between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> column_wise_cross(const matrix_type<3, 3>& m,
                                               const vector3& v) const {
        matrix_type<3, 3> ret;

        vector3 m_col0 = getter::vector<3>(m, 0u, 0u);
        vector3 m_col1 = getter::vector<3>(m, 0u, 1u);
        vector3 m_col2 = getter::vector<3>(m, 0u, 2u);

        getter::set_block(ret, static_cast<vector3>(vector::cross(m_col0, v)),
                          0u, 0u);
        getter::set_block(ret, static_cast<vector3>(vector::cross(m_col1, v)),
                          0u, 1u);
        getter::set_block(ret, static_cast<vector3>(vector::cross(m_col2, v)),
                          0u, 2u);

        return ret;
    }

    /// Column-wise multiplication between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> column_wise_multiply(const matrix_type<3, 3>& m,
                                                  const vector3& v) const {
        matrix_type<3, 3> ret;

        for (std::size_t i = 0u; i < 3u; i++) {
            for (std::size_t j = 0u; j < 3u; j++) {
                getter::element(ret, j, i) = getter::element(m, j, i) * getter::element(v, j);
            }
        }

        return ret;
    }

    /// Cross product matrix
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> cross_matrix(const vector3& v) const {
        matrix_type<3, 3> ret;
        getter::element(ret, 0u, 0u) = 0.f;
        getter::element(ret, 0u, 1u) = -getter::element(v, 2u);
        getter::element(ret, 0u, 2u) = getter::element(v, 1u);
        getter::element(ret, 1u, 0u) = getter::element(v, 2u);
        getter::element(ret, 1u, 1u) = 0.f;
        getter::element(ret, 1u, 2u) = -getter::element(v, 0u);
        getter::element(ret, 2u, 0u) = -getter::element(v, 1u);
        getter::element(ret, 2u, 1u) = getter::element(v, 0u);
        getter::element(ret, 2u, 2u) = 0.f;

        return ret;
    }

    /// Outer product operation
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> outer_product(const vector3& v1,
                                           const vector3& v2) const {
        matrix_type<3, 1> m1;
        getter::element(m1, 0u, 0u) = getter::element(v1, 0u);
        getter::element(m1, 1u, 0u) = getter::element(v1, 1u);
        getter::element(m1, 2u, 0u) = getter::element(v1, 2u);

        matrix_type<1, 3> m2;
        getter::element(m2, 0u, 0u) = getter::element(v2, 0u);
        getter::element(m2, 0u, 1u) = getter::element(v2, 1u);
        getter::element(m2, 0u, 2u) = getter::element(v2, 2u);

        return m1 * m2;
    }

    /// Cholesky decompose
    template <index_type N>
    DETRAY_HOST_DEVICE inline matrix_type<N, N> cholesky_decompose(
        const matrix_type<N, N>& mat) const {

        matrix_type<N, N> L = matrix::zero<matrix_type<N, N>>();

        // Cholesky–Banachiewicz algorithm
        for (std::size_t i = 0u; i < N; i++) {
            for (std::size_t j = 0u; j <= i; j++) {
                scalar_type sum = 0.f;
                for (std::size_t k = 0u; k < j; k++)
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
