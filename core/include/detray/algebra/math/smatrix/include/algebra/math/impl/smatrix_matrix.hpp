/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/smatrix.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

namespace algebra::smatrix::math {

/// Create zero matrix
template <concepts::matrix matrix_t>
ALGEBRA_HOST_DEVICE constexpr matrix_t zero() {
  return matrix_t();
}

/// Create identity matrix
template <concepts::matrix matrix_t>
ALGEBRA_HOST_DEVICE constexpr matrix_t identity() {
  return matrix_t(ROOT::Math::SMatrixIdentity());
}

/// Set input matrix as zero matrix
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr void set_zero(
    ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m) {
  m = zero<ROOT::Math::SMatrix<scalar_t, ROWS, COLS>>();
}

/// Set input matrix as identity matrix
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr void set_identity(
    ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m) {
  m = identity<ROOT::Math::SMatrix<scalar_t, ROWS, COLS>>();
}

/// Create transpose matrix
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> transpose(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m) {
  return ROOT::Math::Transpose(m);
}

/// @returns the determinant of @param m
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr scalar_t determinant(
    const matrix_type<scalar_t, ROWS, COLS> &m) {

  scalar_t det;
  [[maybe_unused]] bool success = m.Det2(det);

  return det;
}

/// @returns the inverse of @param m
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> inverse(
    const matrix_type<scalar_t, ROWS, COLS> &m) {

  int ifail = 0;
  return m.Inverse(ifail);
}

}  // namespace algebra::smatrix::math
