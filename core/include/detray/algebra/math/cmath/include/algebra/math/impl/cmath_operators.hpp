/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// System include(s).
#include <cstddef>
#include <type_traits>

namespace algebra::cmath {

/// @name Operators on 2-element arrays
/// @{

template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar1_t, concepts::scalar scalar2_t>
ALGEBRA_HOST_DEVICE constexpr array_t<scalar1_t, 2> operator*(
    const array_t<scalar1_t, 2> &a, scalar2_t s) {

  return {a[0] * static_cast<scalar1_t>(s), a[1] * static_cast<scalar1_t>(s)};
}

template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar1_t, concepts::scalar scalar2_t>
ALGEBRA_HOST_DEVICE constexpr array_t<scalar1_t, 2> operator*(
    scalar2_t s, const array_t<scalar1_t, 2> &a) {

  return {static_cast<scalar1_t>(s) * a[0], static_cast<scalar1_t>(s) * a[1]};
}

template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr array_t<scalar_t, 2> operator-(
    const array_t<scalar_t, 2> &a, const array_t<scalar_t, 2> &b) {

  return {a[0] - b[0], a[1] - b[1]};
}

template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr array_t<scalar_t, 2> operator+(
    const array_t<scalar_t, 2> &a, const array_t<scalar_t, 2> &b) {

  return {a[0] + b[0], a[1] + b[1]};
}

/// @}

/// @name Operators on 3-element arrays
/// @{

template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar1_t, concepts::scalar scalar2_t>
ALGEBRA_HOST_DEVICE constexpr array_t<scalar1_t, 3> operator*(
    const array_t<scalar1_t, 3> &a, scalar2_t s) {

  return {a[0] * static_cast<scalar1_t>(s), a[1] * static_cast<scalar1_t>(s),
          a[2] * static_cast<scalar1_t>(s)};
}

template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar1_t, concepts::scalar scalar2_t>
ALGEBRA_HOST_DEVICE constexpr array_t<scalar1_t, 3> operator*(
    scalar2_t s, const array_t<scalar1_t, 3> &a) {

  return {static_cast<scalar1_t>(s) * a[0], static_cast<scalar1_t>(s) * a[1],
          static_cast<scalar1_t>(s) * a[2]};
}

template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr array_t<scalar_t, 3> operator-(
    const array_t<scalar_t, 3> &a, const array_t<scalar_t, 3> &b) {

  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr array_t<scalar_t, 3> operator+(
    const array_t<scalar_t, 3> &a, const array_t<scalar_t, 3> &b) {

  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

/// @}

/// @name Operators on matrix
/// @{

template <typename size_type, template <typename, size_type> class array_t,
          concepts::scalar scalar1_t, concepts::scalar scalar2_t,
          size_type ROWS, size_type COLS>
ALGEBRA_HOST_DEVICE constexpr array_t<array_t<scalar1_t, ROWS>, COLS> operator*(
    const array_t<array_t<scalar1_t, ROWS>, COLS> &a, scalar2_t s) {

  array_t<array_t<scalar1_t, ROWS>, COLS> ret;

  for (size_type j = 0; j < COLS; ++j) {
    for (size_type i = 0; i < ROWS; ++i) {
      ret[j][i] = a[j][i] * static_cast<scalar1_t>(s);
    }
  }

  return ret;
}

template <typename size_type, template <typename, size_type> class array_t,
          concepts::scalar scalar1_t, concepts::scalar scalar2_t,
          size_type ROWS, size_type COLS>
ALGEBRA_HOST_DEVICE constexpr array_t<array_t<scalar1_t, ROWS>, COLS> operator*(
    scalar2_t s, const array_t<array_t<scalar1_t, ROWS>, COLS> &a) {

  array_t<array_t<scalar1_t, ROWS>, COLS> ret;

  for (size_type j = 0; j < COLS; ++j) {
    for (size_type i = 0; i < ROWS; ++i) {
      ret[j][i] = a[j][i] * static_cast<scalar1_t>(s);
    }
  }

  return ret;
}

template <typename size_type, template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type M, size_type N, size_type O>
ALGEBRA_HOST_DEVICE constexpr array_t<array_t<scalar_t, M>, O> operator*(
    const array_t<array_t<scalar_t, M>, N> &A,
    const array_t<array_t<scalar_t, N>, O> &B) {

  array_t<array_t<scalar_t, M>, O> C;

  for (size_type j = 0; j < O; ++j) {
    for (size_type i = 0; i < M; ++i) {
      C[j][i] = 0.f;
    }
  }

  for (size_type i = 0; i < N; ++i) {
    for (size_type j = 0; j < O; ++j) {
      for (size_type k = 0; k < M; ++k) {
        C[j][k] += A[i][k] * B[j][i];
      }
    }
  }

  return C;
}

template <typename size_type, template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type ROWS, size_type COLS>
ALGEBRA_HOST_DEVICE constexpr array_t<array_t<scalar_t, ROWS>, COLS> operator+(
    const array_t<array_t<scalar_t, ROWS>, COLS> &A,
    const array_t<array_t<scalar_t, ROWS>, COLS> &B) {

  array_t<array_t<scalar_t, ROWS>, COLS> C;

  for (size_type j = 0; j < COLS; ++j) {
    for (size_type i = 0; i < ROWS; ++i) {
      C[j][i] = A[j][i] + B[j][i];
    }
  }

  return C;
}

template <typename size_type, template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type ROWS, size_type COLS>
ALGEBRA_HOST_DEVICE constexpr array_t<array_t<scalar_t, ROWS>, COLS> operator-(
    const array_t<array_t<scalar_t, ROWS>, COLS> &A,
    const array_t<array_t<scalar_t, ROWS>, COLS> &B) {

  array_t<array_t<scalar_t, ROWS>, COLS> C;

  for (size_type j = 0; j < COLS; ++j) {
    for (size_type i = 0; i < ROWS; ++i) {
      C[j][i] = A[j][i] - B[j][i];
    }
  }

  return C;
}

/// @}

/// @name Operators on matrix * vector
/// @{

template <typename size_type, template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type ROWS, size_type COLS>
ALGEBRA_HOST_DEVICE constexpr array_t<scalar_t, ROWS> operator*(
    const array_t<array_t<scalar_t, ROWS>, COLS> &a,
    const array_t<scalar_t, COLS> &b) {

  array_t<scalar_t, ROWS> ret{0};

  for (size_type j = 0; j < COLS; ++j) {
    for (size_type i = 0; i < ROWS; ++i) {
      ret[i] += a[j][i] * b[j];
    }
  }

  return ret;
}

/// @}

}  // namespace algebra::cmath
