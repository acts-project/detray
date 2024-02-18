/** Detray plugins library, part of the ACTS project
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

namespace detray {

template <typename matrix_operator_t>
struct matrix_helper {

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

    /// Column-wise cross product between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> column_wise_cross(const matrix_type<3, 3>& m,
                                               const vector3& v) const {
        matrix_type<3, 3> ret;

        auto m_col0 = matrix_operator().template block<3, 1>(m, 0, 0);
        auto m_col1 = matrix_operator().template block<3, 1>(m, 0, 1);
        auto m_col2 = matrix_operator().template block<3, 1>(m, 0, 2);

        matrix_operator().set_block(ret, vector::cross(m_col0, v), 0, 0);
        matrix_operator().set_block(ret, vector::cross(m_col1, v), 0, 1);
        matrix_operator().set_block(ret, vector::cross(m_col2, v), 0, 2);

        return ret;
    }

    /// Column-wise multiplication between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> column_wise_multiply(const matrix_type<3, 3>& m,
                                                  const vector3& v) const {
        matrix_type<3, 3> ret;

        for (size_type i = 0; i < 3; i++) {
            for (size_type j = 0; j < 3; j++) {
                matrix_operator().element(ret, j, i) =
                    matrix_operator().element(m, j, i) * v[j];
            }
        }

        return ret;
    }

    /// Cross product matrix
    ///
    /// For a given vector (v) = (x,y,z), return the following matrix:
    ///
    /// [ 0  -z   y ]
    /// [ z   0  -x ]
    /// [-y   x   0 ]
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> cross_matrix(const vector3& v) const {
        matrix_type<3, 3> ret;
        matrix_operator().element(ret, 0, 0) = 0;
        matrix_operator().element(ret, 0, 1) = -v[2];
        matrix_operator().element(ret, 0, 2) = v[1];
        matrix_operator().element(ret, 1, 0) = v[2];
        matrix_operator().element(ret, 1, 1) = 0;
        matrix_operator().element(ret, 1, 2) = -v[0];
        matrix_operator().element(ret, 2, 0) = -v[1];
        matrix_operator().element(ret, 2, 1) = v[0];
        matrix_operator().element(ret, 2, 2) = 0;

        return ret;
    }

    /// Outer product operation
    ///
    /// @return 3x3 matrix M = v * v^T where v is a 3-elements column matrix
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> outer_product(const vector3& v1,
                                           const vector3& v2) const {
        matrix_type<3, 1> m1;
        matrix_operator().element(m1, 0, 0) = v1[0];
        matrix_operator().element(m1, 1, 0) = v1[1];
        matrix_operator().element(m1, 2, 0) = v1[2];

        matrix_type<1, 3> m2;
        matrix_operator().element(m2, 0, 0) = v2[0];
        matrix_operator().element(m2, 0, 1) = v2[1];
        matrix_operator().element(m2, 0, 2) = v2[2];

        return m1 * m2;
    }

    /// Cholesky decompose
    ///
    /// @return the lower triangle matrix (L) that satisfies M = L * L^T
    template <size_type N>
    DETRAY_HOST_DEVICE inline matrix_type<N, N> cholesky_decompose(
        const matrix_type<N, N>& mat) const {

        matrix_type<N, N> L = matrix_operator().template zero<N, N>();

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

    /// @return the rotation matrix (R) that converts the local curvilinear
    /// cartesian to the global cartesian
    ///
    /// The local z-axis (w) is defined as the given vector. If we say that u
    /// and v are the local x and y axis, R is defined as follows:
    ///
    ///     [ u₁  v₁  w₁ ]
    /// R = [ u₂  v₂  w₂ ]
    ///     [ u₃  v₃  w₃ ]
    ///
    /// One way to obtain R is solving Rodrigues' rotation_formula (This is the
    /// formula implemented in axis_rotation.hpp):
    ///
    /// cos_α I + sin_α cross_matrix(k) + (1 - cos_α) * outer_product(k)
    ///
    /// k: normalization of (0,0,1) X w
    /// α: the angle required to rotate (0,0,1) to w around the axis of k
    /// cross_matrix: function defined above
    /// outer_product: function defined_above
    ///
    ///
    /// In a nutshell, the rest of R elements is calculated as:
    ///
    /// u₁ = (w₁ w₁ w₃ + w₂ w₂) / ( 1 - w₃ w₃ )
    /// u₂ = -w₁w₂/(1 + w₃)
    /// u₃ = -w₁
    /// v₁ = -w₁w₂/(1 + w₃)
    /// v₂ = (w₁ w₁ + w₂ w₂ w₃) / ( 1 - w₃ w₃ )
    /// v₃ = -w₂
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> curvilinear_to_global(const vector3& w) const {

        matrix_type<3, 3> R = matrix_operator().template zero<3, 3>();

        const scalar_type w1w2 = w[0] * w[1];

        // Set u
        getter::element(R, 0u, 0u) =
            (w[0] * w[0] * w[2] + w[1] * w[1]) / (1 - w[2] * w[2]);
        getter::element(R, 1u, 0u) = -w1w2 / (1 + w[2]);
        getter::element(R, 2u, 0u) = -w[0];

        // Set v
        getter::element(R, 0u, 1u) = getter::element(R, 1u, 0u);
        getter::element(R, 1u, 1u) =
            (w[0] * w[0] + w[1] * w[1] * w[2]) / (1 - w[2] * w[2]);
        getter::element(R, 2u, 1u) = -w[1];

        // Set w
        getter::element(R, 0u, 2u) = w[0u];
        getter::element(R, 1u, 2u) = w[1u];
        getter::element(R, 2u, 2u) = w[2u];

        return R;
    }

    /// @return the rotation matrix (R) that converts the global cartesian to
    /// the local curvilinear cartesian
    ///
    /// It is the transpose of curvlinear_to_global
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> global_to_curvilinear(const vector3& w) {

        matrix_type<3, 3> R = matrix_operator().template zero<3, 3>();

        const scalar_type w1w2 = w[0] * w[1];

        // Set u
        getter::element(R, 0u, 0u) =
            (w[0] * w[0] * w[2] + w[1] * w[1]) / (1 - w[2] * w[2]);
        getter::element(R, 1u, 0u) = -w1w2 / (1 + w[2]);
        getter::element(R, 2u, 0u) = w[0];

        // Set v
        getter::element(R, 0u, 1u) = getter::element(R, 1u, 0u);
        getter::element(R, 1u, 1u) =
            (w[0] * w[0] + w[1] * w[1] * w[2]) / (1 - w[2] * w[2]);
        getter::element(R, 2u, 1u) = w[1];

        // Set w
        getter::element(R, 0u, 2u) = -w[0u];
        getter::element(R, 1u, 2u) = -w[1u];
        getter::element(R, 2u, 2u) = w[2u];

        return R;
    }
};

}  // namespace detray
