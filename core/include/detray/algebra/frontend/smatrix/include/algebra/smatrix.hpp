/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/generic.hpp"
#include "algebra/math/smatrix.hpp"
#include "algebra/storage/smatrix.hpp"

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::smatrix::matrix_type
/// @{

using smatrix::storage::block;
using smatrix::storage::element;
using smatrix::storage::set_block;
using smatrix::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::cross;
using smatrix::math::dot;
using smatrix::math::eta;
using smatrix::math::norm;
using smatrix::math::normalize;
using smatrix::math::perp;
using smatrix::math::phi;
using smatrix::math::theta;

/// @}

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::determinant;
using smatrix::math::identity;
using smatrix::math::inverse;
using smatrix::math::set_identity;
using smatrix::math::set_zero;
using smatrix::math::transpose;
using smatrix::math::zero;

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

namespace smatrix {

/// @name SMatrix based transforms on @c algebra::smatrix::storage_type
/// @{

template <typename T>
using transform3 = math::transform3<T>;

/// @}

}  // namespace smatrix

namespace plugin {

/// Define the plugin types
/// @{
template <concepts::value V>
struct smatrix {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using size_type = algebra::smatrix::size_type;
  using transform3D = algebra::smatrix::transform3<value_type>;
  using point2D = algebra::smatrix::point2<value_type>;
  using point3D = algebra::smatrix::point3<value_type>;
  using vector3D = algebra::smatrix::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::smatrix::matrix_type<value_type, ROWS, COLS>;
};
/// @}

}  // namespace plugin

}  // namespace algebra
