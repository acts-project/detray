/** Algebra plugins, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"
#include "algebra/storage/fastor.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

namespace algebra::fastor::math {

/// Create zero matrix
template <concepts::matrix matrix_t>
ALGEBRA_HOST_DEVICE constexpr matrix_t zero() {
  return matrix_t(0);
}

/// Create identity matrix
template <concepts::matrix matrix_t>
ALGEBRA_HOST_DEVICE constexpr matrix_t identity() {
  using scalar_t = algebra::traits::value_t<matrix_t>;
  constexpr auto rows{algebra::traits::rows<matrix_t>};
  constexpr auto cols{algebra::traits::columns<matrix_t>};

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
ALGEBRA_HOST_DEVICE constexpr void set_zero(
    matrix_type<scalar_t, ROWS, COLS> &m) {
  m.zeros();
}

/// Set input matrix as identity matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr void set_identity(
    matrix_type<scalar_t, ROWS, COLS> &m) {

  m = identity<matrix_type<scalar_t, ROWS, COLS>>();
}

/// Create transpose matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> transpose(
    const matrix_type<scalar_t, ROWS, COLS> &m) {

  return Fastor::transpose(m);
}

/// @returns the determinant of @param m
template <std::size_t N, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr scalar_t determinant(
    const matrix_type<scalar_t, N, N> &m) {

  return Fastor::determinant(m);
}

/// @returns the inverse of @param m
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> inverse(
    const matrix_type<scalar_t, ROWS, COLS> &m) {

  return Fastor::inverse(m);
}

}  // namespace algebra::fastor::math
