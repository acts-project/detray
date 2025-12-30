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

// ROOT/Smatrix include(s).
#include <Math/Expression.h>
#include <Math/Functions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TMath.h>

// System include(s).
#include <cassert>

namespace algebra::smatrix::storage {

/// Functor used to access elements of SMatrix matrices
struct element_getter {

  template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE constexpr scalar_t &operator()(
      ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, unsigned int row,
      unsigned int col) const {
    assert(row < ROWS);
    assert(col < COLS);
    return m(row, col);
  }

  template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE constexpr scalar_t operator()(
      const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, unsigned int row,
      unsigned int col) const {
    assert(row < ROWS);
    assert(col < COLS);
    return m(row, col);
  }

  template <unsigned int N, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE constexpr scalar_t &operator()(
      ROOT::Math::SMatrix<scalar_t, N, 1> &m, unsigned int row) const {
    assert(row < N);
    return m(row);
  }

  template <unsigned int N, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE constexpr scalar_t operator()(
      const ROOT::Math::SMatrix<scalar_t, N, 1> &m, unsigned int row) const {
    assert(row < N);
    return m(row, 0);
  }

  template <concepts::scalar scalar_t, unsigned int N>
  ALGEBRA_HOST_DEVICE constexpr scalar_t &operator()(
      ROOT::Math::SVector<scalar_t, N> &m, unsigned int row) const {
    assert(row < N);
    return m(row);
  }

  template <unsigned int N, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE constexpr scalar_t operator()(
      const ROOT::Math::SVector<scalar_t, N> &m, unsigned int row) const {
    assert(row < N);
    return m(row);
  }
};  // element_getter

/// Function extracting an element from a matrix (const)
template <concepts::scalar scalar_t, unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE constexpr scalar_t element(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, std::size_t row,
    std::size_t col) {
  return element_getter()(m, static_cast<unsigned int>(row),
                          static_cast<unsigned int>(col));
}

/// Function extracting an element from a matrix (non-const)
template <concepts::scalar scalar_t, unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE constexpr scalar_t &element(
    ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, std::size_t row,
    std::size_t col) {
  return element_getter()(m, static_cast<unsigned int>(row),
                          static_cast<unsigned int>(col));
}

/// Function extracting an element from a matrix (const)
template <concepts::scalar scalar_t, unsigned int N>
ALGEBRA_HOST_DEVICE constexpr scalar_t element(
    const ROOT::Math::SMatrix<scalar_t, N, 1> &m, std::size_t row) {
  return element_getter()(m, static_cast<unsigned int>(row));
}

/// Function extracting an element from a matrix (non-const)
template <concepts::scalar scalar_t, unsigned int N>
ALGEBRA_HOST_DEVICE constexpr scalar_t &element(
    ROOT::Math::SMatrix<scalar_t, N, 1> &m, std::size_t row) {
  return element_getter()(m, static_cast<unsigned int>(row));
}

/// Function extracting an element from a matrix (const)
template <concepts::scalar scalar_t, unsigned int N>
ALGEBRA_HOST_DEVICE constexpr scalar_t element(
    const ROOT::Math::SVector<scalar_t, N> &m, std::size_t row) {
  return element_getter()(m, static_cast<unsigned int>(row));
}

/// Function extracting an element from a matrix (non-const)
template <concepts::scalar scalar_t, unsigned int N>
ALGEBRA_HOST_DEVICE constexpr scalar_t &element(
    ROOT::Math::SVector<scalar_t, N> &m, std::size_t row) {
  return element_getter()(m, static_cast<unsigned int>(row));
}

/// Functor used to extract a block from SMatrix matrices
struct block_getter {

  template <unsigned int ROWS, unsigned int COLS, unsigned int oROWS,
            unsigned int oCOLS, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE ROOT::Math::SMatrix<scalar_t, ROWS, COLS> operator()(
      const ROOT::Math::SMatrix<scalar_t, oROWS, oCOLS> &m, unsigned int row,
      unsigned int col) const {

    return m.template Sub<ROOT::Math::SMatrix<scalar_t, ROWS, COLS>>(row, col);
  }

  template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS,
            concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE ROOT::Math::SVector<scalar_t, SIZE> vector(
      const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, unsigned int row,
      unsigned int col) const {

    // TODO: SMatrix bug?
    // return m.template SubCol<ROOT::Math::SVector<scalar_t, SIZE>>(col, row);

    assert(col < COLS);
    assert(row + SIZE <= ROWS);

    ROOT::Math::SVector<scalar_t, SIZE> ret;

    for (unsigned int irow = row; irow < row + SIZE; ++irow) {
      ret[irow - row] = m(irow, col);
    }

    return ret;
  }
};  // struct block_getter

/// Operator getting a block of a const matrix
template <unsigned int ROWS, unsigned int COLS, unsigned int oROWS,
          unsigned int oCOLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE ROOT::Math::SMatrix<scalar_t, ROWS, COLS> block(
    const ROOT::Math::SMatrix<scalar_t, oROWS, oCOLS> &m, std::size_t row,
    std::size_t col) {
  return block_getter{}.template operator()<ROWS, COLS>(
      m, static_cast<unsigned int>(row), static_cast<unsigned int>(col));
}

/// Function extracting a slice from the matrix used by
/// @c algebra::smatrix::transform3
template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS,
          concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr auto vector(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, std::size_t row,
    std::size_t col) {

  return block_getter{}.template vector<SIZE>(m, static_cast<unsigned int>(row),
                                              static_cast<unsigned int>(col));
}

/// Operator setting a block with a matrix
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t,
          concepts::matrix input_matrix_type>
ALGEBRA_HOST_DEVICE void set_block(
    input_matrix_type &m, const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &b,
    std::size_t row, std::size_t col) {
  for (unsigned int i = 0; i < ROWS; ++i) {
    for (unsigned int j = 0; j < COLS; ++j) {
      m(i + static_cast<unsigned int>(row),
        j + static_cast<unsigned int>(col)) = b(i, j);
    }
  }
}

/// Operator setting a block with a vector
template <unsigned int ROWS, concepts::scalar scalar_t,
          concepts::matrix input_matrix_type>
ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                   const ROOT::Math::SVector<scalar_t, ROWS> &b,
                                   unsigned int row, unsigned int col) {
  for (unsigned int i = 0; i < ROWS; ++i) {
    m(i + row, col) = b[i];
  }
}

}  // namespace algebra::smatrix::storage
