/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/type_traits.hpp"

namespace algebra::generic::matrix::inverse {

/// "Partial Pivot LU Decomposition", assuming a N X N matrix
template <concepts::square_matrix matrix_t, class element_getter_t>
struct partial_pivot_lud {

  using scalar_type = algebra::traits::value_t<matrix_t>;
  using size_type = algebra::traits::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  using decomposition_t =
      typename algebra::generic::matrix::decomposition::partial_pivot_lud<
          matrix_t, element_getter_t>;

  ALGEBRA_HOST_DEVICE constexpr matrix_t operator()(const matrix_t& m) const {

    constexpr size_type N{algebra::traits::rank<matrix_t>};

    const typename decomposition_t::template lud<N> decomp_res =
        decomposition_t()(m);

    // Get the LU decomposition matrix equal to (L - I) + U
    const auto& lu = decomp_res.lu;

    // Permutation vector
    const auto& P = decomp_res.P;

    // Inverse matrix
    matrix_t inv;

    // Calculate inv(A) = inv(U) * inv(L) * P;
    for (size_type j = 0; j < N; j++) {
      for (size_type i = 0; i < N; i++) {
        element_getter_t()(inv, i, j) = static_cast<size_type>(P[i]) == j
                                            ? static_cast<scalar_type>(1.0)
                                            : static_cast<scalar_type>(0.0);

        for (size_type k = 0; k < i; k++) {
          element_getter_t()(inv, i, j) -=
              element_getter_t()(lu, i, k) * element_getter_t()(inv, k, j);
        }
      }

      for (size_type i = N - 1; int(i) >= 0; i--) {
        for (size_type k = i + 1; k < N; k++) {
          element_getter_t()(inv, i, j) -=
              element_getter_t()(lu, i, k) * element_getter_t()(inv, k, j);
        }
        element_getter_t()(inv, i, j) /= element_getter_t()(lu, i, i);
      }
    }

    return inv;
  }
};

}  // namespace algebra::generic::matrix::inverse
