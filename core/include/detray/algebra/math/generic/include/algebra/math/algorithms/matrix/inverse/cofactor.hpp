/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/algorithms/matrix/determinant/cofactor.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace algebra::generic::matrix {

namespace adjoint {

/// "Adjoint getter", assuming a N X N matrix
template <concepts::square_matrix matrix_t, class element_getter_t>
struct cofactor {

  using scalar_type = algebra::traits::value_t<matrix_t>;
  using size_type = algebra::traits::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  ALGEBRA_HOST_DEVICE constexpr matrix_t operator()(const matrix_t &m) const {
    return adjoint_getter_helper<algebra::traits::rank<matrix_t>>()(m);
  }

  template <size_type N, typename Enable = void>
  struct adjoint_getter_helper;

  template <size_type N>
  struct adjoint_getter_helper<N, typename std::enable_if_t<N == 1>> {
    ALGEBRA_HOST_DEVICE constexpr matrix_t operator()(
        const matrix_t & /*m*/) const {
      matrix_t ret;
      element_getter()(ret, 0, 0) = 1;
      return ret;
    }
  };

  template <size_type N>
  struct adjoint_getter_helper<N, typename std::enable_if_t<N != 1>> {

    using determinant_getter =
        determinant::cofactor<matrix_t, element_getter_t>;

    ALGEBRA_HOST_DEVICE constexpr matrix_t operator()(const matrix_t &m) const {

      matrix_t adj;

      // temp is used to store cofactors of m
      int sign = 1;

      // To store cofactors
      matrix_t temp;

      for (size_type i = 0; i < N; i++) {
        for (size_type j = 0; j < N; j++) {
          // Get cofactor of m[i][j]
          typename determinant_getter::template determinant_getter_helper<N>()
              .get_cofactor(m, temp, i, j);

          // sign of adj[j][i] positive if sum of row
          // and column indexes is even.
          sign = ((i + j) % 2 == 0) ? 1 : -1;

          // Interchanging rows and columns to get the
          // transpose of the cofactor matrix
          element_getter()(adj, j, i) =
              sign *
              typename determinant_getter::template determinant_getter_helper<
                  N - 1>()(temp);
        }
      }

      return adj;
    }
  };
};

}  // namespace adjoint

namespace inverse {

/// "inverse getter", assuming a N X N matrix
template <class matrix_t, class element_getter_t>
struct cofactor {

  using scalar_type = algebra::traits::value_t<matrix_t>;
  using size_type = algebra::traits::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  using determinant_getter = determinant::cofactor<matrix_t, element_getter_t>;

  using adjoint_getter = adjoint::cofactor<matrix_t, element_getter_t>;

  ALGEBRA_HOST_DEVICE constexpr matrix_t operator()(const matrix_t &m) const {

    constexpr size_type N{algebra::traits::rank<matrix_t>};

    matrix_t ret;

    // Find determinant of A
    scalar_type det = determinant_getter()(m);

    // TODO: handle singular matrix error
    // if (det == 0) {
    // return ret;
    //}

    auto adj = adjoint_getter()(m);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (size_type i = 0; i < N; i++) {
      for (size_type j = 0; j < N; j++) {
        element_getter()(ret, j, i) = element_getter()(adj, j, i) / det;
      }
    }

    return ret;
  }
};

}  // namespace inverse

}  // namespace algebra::generic::matrix
