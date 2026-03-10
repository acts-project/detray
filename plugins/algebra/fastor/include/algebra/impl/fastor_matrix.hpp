/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/impl/fastor_types.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

namespace detray::algebra::fastor::math {

/// Create zero matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr matrix_t zero() {
    return matrix_t(0);
}

/// Create identity matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr matrix_t identity() {
    using scalar_t = detray::traits::value_t<matrix_t>;
    constexpr auto rows{detray::traits::rows<matrix_t>};
    constexpr auto cols{detray::traits::columns<matrix_t>};

    if constexpr (rows >= cols) {
        matrix_type<scalar_t, rows, rows> identity_matrix;
        identity_matrix.eye2();
        return matrix_t(
            identity_matrix(Fastor::fseq<0, rows>(), Fastor::fseq<0, cols>()));
    } else {
        matrix_type<scalar_t, cols, cols> identity_matrix;
        identity_matrix.eye2();
        return matrix_t(
            identity_matrix(Fastor::fseq<0, rows>(), Fastor::fseq<0, cols>()));
    }
}

/// Set input matrix as zero matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr void set_zero(
    matrix_type<scalar_t, ROWS, COLS> &m) {
    m.zeros();
}

/// Set input matrix as identity matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr void set_identity(
    matrix_type<scalar_t, ROWS, COLS> &m) {

    m = identity<matrix_type<scalar_t, ROWS, COLS>>();
}

/// Create transpose matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> transpose(
    const matrix_type<scalar_t, ROWS, COLS> &m) {

    return Fastor::transpose(m);
}

/// Column-wise cross product
/// @{
template <std::size_t ROWS, std::size_t COLS, std::size_t N,
          concepts::scalar scalar_t, typename derived_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS>
column_wise_cross(const matrix_type<scalar_t, COLS, ROWS> &m,
                  const Fastor::AbstractTensor<derived_t, N> &v) {
    return m * v;
}

template <std::size_t N, std::size_t M, typename derived_1_t,
          typename derived_2_t>
DETRAY_HOST_DEVICE constexpr auto column_wise_cross(
    const Fastor::AbstractTensor<derived_1_t, M> &m,
    const Fastor::AbstractTensor<derived_2_t, N> &v) {
    return m * v;
}
/// @}

/// Column-wise product
/// @{
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS>
column_wise_multiply(const matrix_type<scalar_t, COLS, ROWS> &m,
                     const matrix_type<scalar_t, COLS, ROWS> &b) {
    return m * b;
}

template <std::size_t N, std::size_t M, typename derived_1_t,
          typename derived_2_t>
DETRAY_HOST_DEVICE constexpr auto column_wise_multiply(
    const Fastor::AbstractTensor<derived_1_t, M> &m,
    const Fastor::AbstractTensor<derived_2_t, N> &b) {
    return m * b;
}
/// @}

/// @returns the determinant of @param m
template <std::size_t N, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr scalar_t determinant(
    const matrix_type<scalar_t, N, N> &m) {

    return Fastor::determinant(m);
}

/// @returns the inverse of @param m
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> inverse(
    const matrix_type<scalar_t, ROWS, COLS> &m) {

    return Fastor::inverse(m);
}

}  // namespace detray::algebra::fastor::math
