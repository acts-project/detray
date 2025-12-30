/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

// System include(s).
#include <array>

namespace algebra::generic::matrix::decomposition {

/// "Partial Pivot LU Decomposition", assuming a N X N matrix
template <concepts::matrix matrix_t, class element_getter_t>
struct partial_pivot_lud {

  using scalar_type = algebra::traits::value_t<matrix_t>;
  using size_type = algebra::traits::index_t<matrix_t>;
  using vector_type = algebra::traits::vector_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  template <size_type N>
  struct lud {
    // LU decomposition matrix, equal to (L - I) + U, where the diagonal
    // components of L is always 1
    matrix_t lu;

    // Permutation vector
    vector_type P;

    // Number of pivots
    int n_pivot = 0;
  };

  ALGEBRA_HOST_DEVICE constexpr lud<algebra::traits::rank<matrix_t>> operator()(
      const matrix_t& m) const {

    constexpr size_type N{algebra::traits::rank<matrix_t>};

    // LU decomposition matrix
    matrix_t lu = m;

    // Permutation
    vector_type P;

    // Max index and value
    size_type max_idx;
    scalar_type max_val;
    scalar_type abs_val;

    // Number of pivoting
    int n_pivot = N;

    // Rows for swapping
    vector_type row_0;
    vector_type row_1;

    // Unit permutation matrix, P[N] initialized with N
    for (size_type i = 0; i < N; i++) {
      P[i] = static_cast<scalar_type>(i);
    }

    for (size_type i = 0; i < N; i++) {
      max_val = 0;
      max_idx = i;

      for (size_type k = i; k < N; k++) {
        abs_val = algebra::math::fabs(element_getter()(lu, k, i));

        if (abs_val > max_val) {

          max_val = abs_val;
          max_idx = k;
        }
      }

      if (max_idx != i) {
        // Pivoting P
        auto j = P[i];

        P[i] = P[max_idx];
        P[max_idx] = j;

        // Pivoting rows of A
        for (size_type q = 0; q < N; q++) {
          row_0[q] = element_getter_t()(lu, i, q);
          row_1[q] = element_getter_t()(lu, max_idx, q);
        }
        for (size_type q = 0; q < N; q++) {
          element_getter_t()(lu, i, q) = row_1[q];
          element_getter_t()(lu, max_idx, q) = row_0[q];
        }

        // counting pivots starting from N (for determinant)
        n_pivot++;
      }

      for (size_type j = i + 1; j < N; j++) {
        // m[j][i] /= m[i][i];
        element_getter_t()(lu, j, i) /= element_getter_t()(lu, i, i);

        for (size_type k = i + 1; k < N; k++) {
          // m[j][k] -= m[j][i] * m[i][k];
          element_getter_t()(lu, j, k) -=
              element_getter_t()(lu, j, i) * element_getter_t()(lu, i, k);
        }
      }
    }

    return {lu, P, n_pivot};
  }
};

}  // namespace algebra::generic::matrix::decomposition
