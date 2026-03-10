/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/impl/eigen_types.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 20012
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__

namespace detray::algebra::eigen::math {

/// Create zero matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr matrix_t zero() {
    return matrix_t::Zero();
}

/// Create identity matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr matrix_t identity() {
    return matrix_t::Identity();
}

/// Set input matrix as zero matrix
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr void set_zero(Eigen::MatrixBase<derived_type> &m) {
    m.setZero();
}

/// Set input matrix as identity matrix
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr void set_identity(
    Eigen::MatrixBase<derived_type> &m) {
    m.setIdentity();
}

/// Create transpose matrix
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr matrix_type<
    typename Eigen::MatrixBase<derived_type>::value_type,
    Eigen::MatrixBase<derived_type>::ColsAtCompileTime,
    Eigen::MatrixBase<derived_type>::RowsAtCompileTime>
transpose(const Eigen::MatrixBase<derived_type> &m) {
    return m.transpose();
}

/// Column-wise cross product
template <typename derived_type_1, typename derived_type_2>
DETRAY_HOST_DEVICE constexpr matrix_type<
    typename Eigen::MatrixBase<derived_type_1>::value_type,
    Eigen::MatrixBase<derived_type_1>::RowsAtCompileTime,
    Eigen::MatrixBase<derived_type_1>::ColsAtCompileTime>
column_wise_cross(const Eigen::MatrixBase<derived_type_1> &m,
                  const Eigen::MatrixBase<derived_type_2> &v) {
    return m.colwise().cross(v);
}

/// Column-wise product
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr matrix_type<
    typename Eigen::MatrixBase<derived_type>::value_type,
    Eigen::MatrixBase<derived_type>::RowsAtCompileTime,
    Eigen::MatrixBase<derived_type>::ColsAtCompileTime>
column_wise_multiply(const Eigen::MatrixBase<derived_type> &m,
                     const Eigen::MatrixBase<derived_type> &b) {
    return m.cwiseProduct(b);
}

/// @returns the determinant of @param m
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr
    typename Eigen::MatrixBase<derived_type>::value_type
    determinant(const Eigen::MatrixBase<derived_type> &m) {
    return m.determinant();
}

/// @returns the inverse of @param m
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr matrix_type<
    typename Eigen::MatrixBase<derived_type>::value_type,
    Eigen::MatrixBase<derived_type>::RowsAtCompileTime,
    Eigen::MatrixBase<derived_type>::ColsAtCompileTime>
inverse(const Eigen::MatrixBase<derived_type> &m) {
    return m.inverse();
}

/// @returns the fatcor L in the decomposition of @param m
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr matrix_type<
    typename Eigen::MatrixBase<derived_type>::value_type,
    Eigen::MatrixBase<derived_type>::RowsAtCompileTime,
    Eigen::MatrixBase<derived_type>::ColsAtCompileTime>
cholesky_decomposition(const Eigen::MatrixBase<derived_type> &m) {
    Eigen::LLT<Eigen::MatrixBase<derived_type>> LL_T(m);
    return LL_T.matrixL();
}

}  // namespace detray::algebra::eigen::math
