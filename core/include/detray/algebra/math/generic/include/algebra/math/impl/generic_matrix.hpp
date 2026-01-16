/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/algorithms/utils/algorithm_finder.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::generic::math {

/// Create zero matrix - generic transform3
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE constexpr M zero() {

  using index_t = algebra::traits::index_t<M>;
  using element_getter_t = algebra::traits::element_getter_t<M>;

  M ret;

  for (index_t j = 0; j < algebra::traits::columns<M>; ++j) {
    for (index_t i = 0; i < algebra::traits::rows<M>; ++i) {
      element_getter_t{}(ret, i, j) = 0;
    }
  }

  return ret;
}

/// Create identity matrix - generic transform3
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE constexpr M identity() {

  using index_t = algebra::traits::index_t<M>;
  using element_getter_t = algebra::traits::element_getter_t<M>;

  auto ret{zero<M>()};

  for (index_t i = 0; i < algebra::traits::rank<M>; ++i) {
    element_getter_t{}(ret, i, i) = 1;
  }

  return ret;
}

/// Set @param m as zero matrix
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE constexpr void set_zero(M &m) {
  m = zero<M>();
}

/// Set @param m as identity matrix
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE constexpr void set_identity(M &m) {
  m = identity<M>();
}

/// @returns the transpose matrix of @param m
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE constexpr auto transpose(const M &m) {

  using index_t = algebra::traits::index_t<M>;
  using value_t = algebra::traits::value_t<M>;
  using element_getter_t = algebra::traits::element_getter_t<M>;

  constexpr index_t rows{algebra::traits::rows<M>};
  constexpr index_t columns{algebra::traits::columns<M>};

  algebra::traits::get_matrix_t<M, columns, rows, value_t> ret;

  for (index_t i = 0; i < rows; ++i) {
    for (index_t j = 0; j < columns; ++j) {
      element_getter_t{}(ret, j, i) = element_getter_t{}(m, i, j);
    }
  }

  return ret;
}

// Set matrix C to the product AB
template <concepts::matrix MC, concepts::matrix MA, concepts::matrix MB>
ALGEBRA_HOST_DEVICE constexpr void
set_product(MC &C, const MA &A, const MB &B) requires(
    algebra::concepts::matrix_multipliable_into<MA, MB, MC>) {
  using index_t = algebra::traits::index_t<MC>;
  using value_t = algebra::traits::value_t<MC>;

  for (index_t i = 0; i < algebra::traits::rows<MC>; ++i) {
    for (index_t j = 0; j < algebra::traits::columns<MC>; ++j) {
      value_t t = 0.f;

      for (index_t k = 0; k < algebra::traits::rows<MB>; ++k) {
        t += algebra::traits::element_getter_t<MA>()(A, i, k) *
             algebra::traits::element_getter_t<MB>()(B, k, j);
      }

      algebra::traits::element_getter_t<MC>()(C, i, j) = t;
    }
  }
}

// Set matrix C to the product A^TB
template <concepts::matrix MC, concepts::matrix MA, concepts::matrix MB>
ALGEBRA_HOST_DEVICE constexpr void
set_product_left_transpose(MC &C, const MA &A, const MB &B) requires(
    algebra::concepts::matrix_multipliable_into<
        decltype(transpose(std::declval<MA>())), MB, MC>) {
  using index_t = algebra::traits::index_t<MC>;
  using value_t = algebra::traits::value_t<MC>;

  for (index_t i = 0; i < algebra::traits::rows<MC>; ++i) {
    for (index_t j = 0; j < algebra::traits::columns<MC>; ++j) {
      value_t t = 0.f;

      for (index_t k = 0; k < algebra::traits::rows<MB>; ++k) {
        t += algebra::traits::element_getter_t<MA>()(A, k, i) *
             algebra::traits::element_getter_t<MB>()(B, k, j);
      }

      algebra::traits::element_getter_t<MC>()(C, i, j) = t;
    }
  }
}

// Set matrix C to the product AB^T
template <concepts::matrix MC, concepts::matrix MA, concepts::matrix MB>
ALGEBRA_HOST_DEVICE constexpr void
set_product_right_transpose(MC &C, const MA &A, const MB &B) requires(
    algebra::concepts::matrix_multipliable_into<
        MA, decltype(transpose(std::declval<MB>())), MC>) {
  using index_t = algebra::traits::index_t<MC>;
  using value_t = algebra::traits::value_t<MC>;

  for (index_t i = 0; i < algebra::traits::rows<MC>; ++i) {
    for (index_t j = 0; j < algebra::traits::columns<MC>; ++j) {
      value_t t = 0.f;

      for (index_t k = 0; k < algebra::traits::columns<MA>; ++k) {
        t += algebra::traits::element_getter_t<MA>()(A, i, k) *
             algebra::traits::element_getter_t<MB>()(B, j, k);
      }

      algebra::traits::element_getter_t<MC>()(C, i, j) = t;
    }
  }
}

// Set matrix A to the product AB in place
template <concepts::matrix MA, concepts::matrix MB>
ALGEBRA_HOST_DEVICE constexpr void
set_inplace_product_right(MA &A, const MB &B) requires(
    algebra::concepts::matrix_multipliable_into<MA, MB, MA>) {
  using index_t = algebra::traits::index_t<MA>;
  using value_t = algebra::traits::value_t<MA>;

  for (index_t i = 0; i < algebra::traits::rows<MA>; ++i) {
    algebra::traits::get_matrix_t<MA, 1, algebra::traits::columns<MA>, value_t>
        Q;

    for (index_t j = 0; j < algebra::traits::columns<MA>; ++j) {
      algebra::traits::element_getter_t<decltype(Q)>()(Q, 0, j) =
          algebra::traits::element_getter_t<MA>()(A, i, j);
    }

    for (index_t j = 0; j < algebra::traits::columns<MA>; ++j) {
      value_t t = 0.f;

      for (index_t k = 0; k < algebra::traits::rows<MB>; ++k) {
        t += algebra::traits::element_getter_t<decltype(Q)>()(Q, 0, k) *
             algebra::traits::element_getter_t<MB>()(B, k, j);
      }

      algebra::traits::element_getter_t<MA>()(A, i, j) = t;
    }
  }
}

// Set matrix A to the product BA in place
template <concepts::matrix MA, concepts::matrix MB>
ALGEBRA_HOST_DEVICE constexpr void
set_inplace_product_left(MA &A, const MB &B) requires(
    algebra::concepts::matrix_multipliable_into<MB, MA, MA>) {
  using index_t = algebra::traits::index_t<MA>;
  using value_t = algebra::traits::value_t<MA>;

  for (index_t j = 0; j < algebra::traits::columns<MA>; ++j) {
    algebra::traits::get_matrix_t<MA, 1, algebra::traits::columns<MA>, value_t>
        Q;

    for (index_t i = 0; i < algebra::traits::rows<MA>; ++i) {
      algebra::traits::element_getter_t<decltype(Q)>()(Q, 0, i) =
          algebra::traits::element_getter_t<MA>()(A, i, j);
    }

    for (index_t i = 0; i < algebra::traits::rows<MA>; ++i) {
      value_t t = 0.f;

      for (index_t k = 0; k < algebra::traits::columns<MB>; ++k) {
        t += algebra::traits::element_getter_t<MB>()(B, i, k) *
             algebra::traits::element_getter_t<decltype(Q)>()(Q, 0, k);
      }

      algebra::traits::element_getter_t<MA>()(A, i, j) = t;
    }
  }
}

// Set matrix A to the product AB^T in place
template <concepts::matrix MA, concepts::matrix MB>
ALGEBRA_HOST_DEVICE constexpr void
set_inplace_product_right_transpose(MA &A, const MB &B) requires(
    algebra::concepts::matrix_multipliable_into<
        MA, decltype(transpose(std::declval<MB>())), MA>) {
  using index_t = algebra::traits::index_t<MA>;
  using value_t = algebra::traits::value_t<MA>;

  for (index_t i = 0; i < algebra::traits::rows<MA>; ++i) {
    algebra::traits::get_matrix_t<MA, 1, algebra::traits::columns<MA>, value_t>
        Q;

    for (index_t j = 0; j < algebra::traits::columns<MA>; ++j) {
      algebra::traits::element_getter_t<decltype(Q)>()(Q, 0, j) =
          algebra::traits::element_getter_t<MA>()(A, i, j);
    }

    for (index_t j = 0; j < algebra::traits::columns<MA>; ++j) {
      value_t T = 0.f;

      for (index_t k = 0; k < algebra::traits::columns<MB>; ++k) {
        T += algebra::traits::element_getter_t<decltype(Q)>()(Q, 0, k) *
             algebra::traits::element_getter_t<MB>()(B, j, k);
      }

      algebra::traits::element_getter_t<MA>()(A, i, j) = T;
    }
  }
}

// Set matrix A to the product B^TA in place
template <concepts::matrix MA, concepts::matrix MB>
ALGEBRA_HOST_DEVICE constexpr void
set_inplace_product_left_transpose(MA &A, const MB &B) requires(
    algebra::concepts::matrix_multipliable_into<
        decltype(transpose(std::declval<MB>())), MA, MA>) {
  using index_t = algebra::traits::index_t<MA>;
  using value_t = algebra::traits::value_t<MA>;

  for (index_t j = 0; j < algebra::traits::columns<MA>; ++j) {
    algebra::traits::get_matrix_t<MA, 1, algebra::traits::columns<MA>, value_t>
        Q;

    for (index_t i = 0; i < algebra::traits::rows<MA>; ++i) {
      algebra::traits::element_getter_t<decltype(Q)>()(Q, 0, i) =
          algebra::traits::element_getter_t<MA>()(A, i, j);
    }

    for (index_t i = 0; i < algebra::traits::rows<MA>; ++i) {
      value_t T = 0.f;

      for (index_t k = 0; k < algebra::traits::rows<MB>; ++k) {
        T += algebra::traits::element_getter_t<MB>()(B, k, i) *
             algebra::traits::element_getter_t<decltype(Q)>()(Q, 0, k);
      }

      algebra::traits::element_getter_t<MA>()(A, i, j) = T;
    }
  }
}

template <bool transpose, concepts::matrix M>
ALGEBRA_HOST_DEVICE constexpr algebra::traits::value_t<M> transposable_get(
    const M &A, const algebra::traits::index_t<M> i,
    const algebra::traits::index_t<M> j) {
  if constexpr (transpose) {
    return algebra::traits::element_getter_t<M>()(A, j, i);
  } else {
    return algebra::traits::element_getter_t<M>()(A, i, j);
  }
}

template <bool transpose_left, bool transpose_right, concepts::matrix MA,
          concepts::matrix MB>
ALGEBRA_HOST_DEVICE constexpr auto
transposed_product(const MA &A, const MB &B) requires(
    algebra::concepts::matrix_multipliable<
        std::conditional_t<transpose_left,
                           decltype(transpose(std::declval<MA>())), MA>,
        std::conditional_t<transpose_right,
                           decltype(transpose(std::declval<MB>())), MB> >) {
  using index_t = algebra::traits::index_t<MA>;
  using value_t = algebra::traits::value_t<MA>;

  using effective_lhs_t =
      std::conditional_t<transpose_left,
                         decltype(transpose(std::declval<MA>())), MA>;
  using effective_rhs_t =
      std::conditional_t<transpose_right,
                         decltype(transpose(std::declval<MB>())), MB>;

  using result_t =
      algebra::traits::get_matrix_t<MA, algebra::traits::rows<effective_lhs_t>,
                                    algebra::traits::columns<effective_rhs_t>,
                                    value_t>;

  result_t C = zero<result_t>();

  for (index_t i = 0; i < algebra::traits::columns<effective_lhs_t>; ++i) {
    for (index_t j = 0; j < algebra::traits::columns<result_t>; ++j) {
      for (index_t k = 0; k < algebra::traits::rows<result_t>; ++k) {
        algebra::traits::element_getter_t<result_t>()(C, k, j) +=
            transposable_get<transpose_left>(A, k, i) *
            transposable_get<transpose_right>(B, i, j);
      }
    }
  }

  return C;
}

/// @returns the determinant of @param m
template <concepts::square_matrix M>
ALGEBRA_HOST_DEVICE constexpr algebra::traits::scalar_t<M> determinant(
    const M &m) {

  return determinant_t<M>{}(m);
}

/// @returns the determinant of @param m
template <concepts::square_matrix M>
ALGEBRA_HOST_DEVICE constexpr M inverse(const M &m) {

  return inversion_t<M>{}(m);
}

}  // namespace algebra::generic::math
