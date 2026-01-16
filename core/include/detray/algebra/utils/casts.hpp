/** Algebra plugins, part of the ACTS project
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/concepts.hpp"
#include "algebra/qualifiers.hpp"

// System include(s)
#include <concepts>

namespace detray::algebra {

/// Cast a salar (might be simd) @param s to the precision given by
/// @tparam value_t
template <concepts::value value_t, concepts::scalar scalar_t>
requires std::convertible_to<scalar_t, value_t>
    ALGEBRA_HOST_DEVICE constexpr auto cast_to(const scalar_t& s) {
  return static_cast<value_t>(s);
}

/// Cast a generic vector or point @param v to the precision given by
/// @tparam value_t
template <concepts::value value_t, typename vector_t>
requires(concepts::vector<vector_t> ||
         concepts::point<vector_t>) ALGEBRA_HOST_DEVICE
    constexpr auto cast_to(const vector_t& v) {

  using index_t = algebra::traits::index_t<vector_t>;

  constexpr index_t size{algebra::traits::size<vector_t>};

  using new_vector_t = algebra::traits::get_vector_t<vector_t, size, value_t>;
  new_vector_t ret;

  static_assert(std::same_as<value_t, algebra::traits::value_t<new_vector_t>>);

  for (index_t i = 0; i < size; ++i) {
    ret[i] = ::algebra::cast_to<value_t>(v[i]);
  }

  return ret;
}

/// Cast a column matrix @param v to the precision given by @tparam value_t
template <concepts::value value_t, concepts::column_matrix vector_t>
ALGEBRA_HOST_DEVICE constexpr auto cast_to(const vector_t& v) {

  using index_t = algebra::traits::index_t<vector_t>;
  using element_getter_t = algebra::traits::element_getter_t<vector_t>;

  constexpr index_t rows{algebra::traits::rows<vector_t>};

  using new_vector_t =
      algebra::traits::get_matrix_t<vector_t, rows, 1, value_t>;
  new_vector_t ret;

  static_assert(std::same_as<value_t, algebra::traits::value_t<new_vector_t>>);

  for (index_t i = 0; i < rows; ++i) {
    element_getter_t{}(ret, i, 0) =
        ::algebra::cast_to<value_t>(element_getter_t{}(v, i, 0));
  }

  return ret;
}

/// Cast a generic matrix @param v to the precision given by @tparam value_t
template <concepts::value value_t, concepts::matrix matrix_t>
ALGEBRA_HOST_DEVICE constexpr auto cast_to(const matrix_t& m) {

  using index_t = algebra::traits::index_t<matrix_t>;
  using element_getter_t = algebra::traits::element_getter_t<matrix_t>;

  constexpr index_t rows{algebra::traits::rows<matrix_t>};
  constexpr index_t columns{algebra::traits::columns<matrix_t>};

  using new_matrix_t =
      algebra::traits::get_matrix_t<matrix_t, rows, columns, value_t>;
  new_matrix_t ret;

  static_assert(std::same_as<value_t, algebra::traits::value_t<new_matrix_t>>);

  for (index_t j = 0; j < columns; ++j) {
    for (index_t i = 0; i < rows; ++i) {
      element_getter_t{}(ret, i, j) =
          ::algebra::cast_to<value_t>(element_getter_t{}(m, i, j));
    }
  }

  return ret;
}

/// Cast a 3D transform @param trf to the precision given by @tparam scalar_t
template <concepts::scalar scalar_t, concepts::transform3D transform_t>
requires(!concepts::simd_scalar<scalar_t> &&
         !concepts::simd_scalar<typename transform_t::scalar_type>)
    ALGEBRA_HOST_DEVICE constexpr auto cast_to(const transform_t& trf) {
  using new_trf3_t = typename transform_t::template other_type<scalar_t>;

  return new_trf3_t{::algebra::cast_to<scalar_t>(trf.matrix()),
                    ::algebra::cast_to<scalar_t>(trf.matrix_inverse())};
}

}  // namespace algebra
