/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/algorithms/matrix/determinant/hard_coded.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace algebra::generic::matrix::inverse {

/// "inverse getter", assuming a N X N matrix
template <concepts::square_matrix matrix_t, class element_getter_t>
struct hard_coded {

  using scalar_type = algebra::traits::value_t<matrix_t>;
  using size_type = algebra::traits::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  using determinant_getter =
      determinant::hard_coded<matrix_t, element_getter_t>;

  // 2 X 2 matrix inverse
  template <typename M = matrix_t>
  requires(algebra::traits::rank<M> == 2) ALGEBRA_HOST_DEVICE constexpr matrix_t
  operator()(const matrix_t &m) const {

    matrix_t ret;

    scalar_type det = determinant_getter()(m);

    element_getter()(ret, 0, 0) = element_getter()(m, 1, 1) / det;
    element_getter()(ret, 0, 1) = -1 * element_getter()(m, 0, 1) / det;
    element_getter()(ret, 1, 0) = -1 * element_getter()(m, 1, 0) / det;
    element_getter()(ret, 1, 1) = element_getter()(m, 0, 0) / det;

    return ret;
  }

  // 4 X 4 matrix inverse
  template <typename M = matrix_t>
  requires(algebra::traits::rank<M> == 4) ALGEBRA_HOST_DEVICE constexpr matrix_t
  operator()(const matrix_t &m) const {

    matrix_t ret;
    element_getter()(ret, 0, 0) =
        element_getter()(m, 1, 2) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 1, 3) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 1, 3) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 1, 1) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 1, 2) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 1, 1) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 3);
    element_getter()(ret, 0, 1) =
        element_getter()(m, 0, 3) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 2) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 3) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 1) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 2) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 3) -
        element_getter()(m, 0, 1) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 3);
    element_getter()(ret, 0, 2) =
        element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 3);
    element_getter()(ret, 0, 3) =
        element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 1) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 1) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 2) +
        element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 2) +
        element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 3) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 3);
    element_getter()(ret, 1, 0) =
        element_getter()(m, 1, 3) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 2) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 3) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 1, 0) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 1, 2) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 3) -
        element_getter()(m, 1, 0) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 3);
    element_getter()(ret, 1, 1) =
        element_getter()(m, 0, 2) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 0) +
        element_getter()(m, 0, 3) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 0) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 2) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 0, 0) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 3);
    element_getter()(ret, 1, 2) =
        element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 3) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 3);
    element_getter()(ret, 1, 3) =
        element_getter()(m, 0, 2) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 0) +
        element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 2) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 2) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 3) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 3);
    element_getter()(ret, 2, 0) =
        element_getter()(m, 1, 1) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 3) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 0) +
        element_getter()(m, 1, 3) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 1, 0) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 1, 1) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 1, 0) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 3);
    element_getter()(ret, 2, 1) =
        element_getter()(m, 0, 3) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 1) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 0) * element_getter()(m, 2, 3) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 1) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 3) -
        element_getter()(m, 0, 0) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 3);
    element_getter()(ret, 2, 2) =
        element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 0) +
        element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 3) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 3);
    element_getter()(ret, 2, 3) =
        element_getter()(m, 0, 3) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 0) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 0) -
        element_getter()(m, 0, 3) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 1) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 3) *
            element_getter()(m, 2, 1) +
        element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 3) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 3);
    element_getter()(ret, 3, 0) =
        element_getter()(m, 1, 2) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 1) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 1, 2) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 1, 0) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 1, 1) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 1, 0) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 2);
    element_getter()(ret, 3, 1) =
        element_getter()(m, 0, 1) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 2) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 0) +
        element_getter()(m, 0, 2) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 0) * element_getter()(m, 2, 2) *
            element_getter()(m, 3, 1) -
        element_getter()(m, 0, 1) * element_getter()(m, 2, 0) *
            element_getter()(m, 3, 2) +
        element_getter()(m, 0, 0) * element_getter()(m, 2, 1) *
            element_getter()(m, 3, 2);
    element_getter()(ret, 3, 2) =
        element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 0) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
            element_getter()(m, 3, 1) +
        element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
            element_getter()(m, 3, 2) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
            element_getter()(m, 3, 2);
    element_getter()(ret, 3, 3) =
        element_getter()(m, 0, 1) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 0) -
        element_getter()(m, 0, 2) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 0) +
        element_getter()(m, 0, 2) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 1) -
        element_getter()(m, 0, 0) * element_getter()(m, 1, 2) *
            element_getter()(m, 2, 1) -
        element_getter()(m, 0, 1) * element_getter()(m, 1, 0) *
            element_getter()(m, 2, 2) +
        element_getter()(m, 0, 0) * element_getter()(m, 1, 1) *
            element_getter()(m, 2, 2);

    scalar_type idet = static_cast<scalar_type>(1.) / determinant_getter()(ret);
    for (unsigned int c = 0; c < 4; ++c) {
      for (unsigned int r = 0; r < 4; ++r) {
        element_getter()(ret, c, r) *= idet;
      }
    }
    return ret;
  }
};

}  // namespace algebra::generic::matrix::inverse
