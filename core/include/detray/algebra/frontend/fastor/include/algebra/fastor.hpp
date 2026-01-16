/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/fastor.hpp"
#include "algebra/math/generic.hpp"
#include "algebra/storage/fastor.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::fastor::matrix_type
/// @{

using fastor::storage::block;
using fastor::storage::element;
using fastor::storage::set_block;
using fastor::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::fastor::storage_type
/// @{

using fastor::math::cross;
using fastor::math::dot;
using fastor::math::eta;
using fastor::math::norm;
using fastor::math::normalize;
using fastor::math::perp;
using fastor::math::phi;
using fastor::math::theta;

/// @}

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::fastor::storage_type
/// @{

using fastor::math::determinant;
using fastor::math::identity;
using fastor::math::inverse;
using fastor::math::set_identity;
using fastor::math::set_zero;
using fastor::math::transpose;
using fastor::math::zero;

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

namespace fastor {

/// @name Transform on @c algebra::fastor::storage_type
/// @{

template <concepts::scalar T>
using transform3 = math::transform3<T>;

/// @}

}  // namespace fastor

namespace plugin {

/// Define the plugin types
/// @{
template <concepts::value V>
struct fastor {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using size_type = algebra::fastor::size_type;
  using transform3D = algebra::fastor::transform3<value_type>;
  using point2D = algebra::fastor::point2<value_type>;
  using point3D = algebra::fastor::point3<value_type>;
  using vector3D = algebra::fastor::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::fastor::matrix_type<value_type, ROWS, COLS>;
};
/// @}

}  // namespace plugin

}  // namespace algebra
