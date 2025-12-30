/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/generic.hpp"
#include "algebra/math/vc_aos.hpp"
#include "algebra/storage/vc_aos.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::vc_aos::matrix_type
/// @{

using vc_aos::storage::block;
using vc_aos::storage::element;
using vc_aos::storage::set_block;
using vc_aos::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc_aos types
/// @{

// Vc array specific
using vc_aos::math::cross;
using vc_aos::math::dot;
using vc_aos::math::eta;
using vc_aos::math::norm;
using vc_aos::math::normalize;
using vc_aos::math::perp;
using vc_aos::math::phi;
using vc_aos::math::theta;

/// @}

}  // namespace vector

// Use special algorithms for 4 dimensional matrices
namespace generic {

// Determinant algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct determinant_selector<4, vc_aos::matrix_type<T, ROWS, COLS>,
                            vc_aos::element_getter> {
  using type =
      matrix::determinant::hard_coded<vc_aos::matrix_type<T, ROWS, COLS>,
                                      vc_aos::element_getter>;
};

// Inversion algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct inversion_selector<4, vc_aos::matrix_type<T, ROWS, COLS>,
                          vc_aos::element_getter> {
  using type = matrix::inverse::hard_coded<vc_aos::matrix_type<T, ROWS, COLS>,
                                           vc_aos::element_getter>;
};

}  // namespace generic

namespace matrix {

/// @name Matrix functions on @c algebra::vc_aos types
/// @{

using vc_aos::math::determinant;
using vc_aos::math::identity;
using vc_aos::math::inverse;
using vc_aos::math::set_identity;
using vc_aos::math::set_zero;
using vc_aos::math::transpose;
using vc_aos::math::zero;

using generic::math::set_inplace_product_left;
using generic::math::set_inplace_product_left_transpose;
using generic::math::set_inplace_product_right;
using generic::math::set_inplace_product_right_transpose;
using generic::math::set_product;
using generic::math::set_product_left_transpose;
using generic::math::set_product_right_transpose;
using generic::math::transposed_product;

/// @}

}  // namespace matrix

namespace vc_aos {

/// @name Vc based transforms on @c algebra::vc_aos::storage_type
/// @{

template <concepts::value T>
using transform3 = math::transform3<vc_aos::storage_type, T>;

/// @}

}  // namespace vc_aos

namespace plugin {

/// Define the plugin types
/// @{
template <concepts::value V>
struct vc_aos {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using size_type = algebra::vc_aos::size_type;
  using transform3D = algebra::vc_aos::transform3<value_type>;
  using point2D = algebra::vc_aos::point2<value_type>;
  using point3D = algebra::vc_aos::point3<value_type>;
  using vector3D = algebra::vc_aos::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::vc_aos::matrix_type<value_type, ROWS, COLS>;
};
/// @}

}  // namespace plugin

}  // namespace algebra
