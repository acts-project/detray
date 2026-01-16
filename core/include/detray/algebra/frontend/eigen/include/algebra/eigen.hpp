/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/eigen.hpp"
#include "algebra/math/generic.hpp"
#include "algebra/storage/eigen.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::eigen
/// @{

using eigen::storage::block;
using eigen::storage::element;
using eigen::storage::set_block;
using eigen::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::eigen::storage_type
/// @{

using eigen::math::cross;
using eigen::math::dot;
using eigen::math::eta;
using eigen::math::norm;
using eigen::math::normalize;
using eigen::math::perp;
using eigen::math::phi;
using eigen::math::theta;

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::eigen::storage_type
/// @{

using eigen::math::determinant;
using eigen::math::identity;
using eigen::math::inverse;
using eigen::math::set_identity;
using eigen::math::set_zero;
using eigen::math::transpose;
using eigen::math::zero;

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

namespace eigen {

template <concepts::scalar T>
using transform3 = math::transform3<T>;

}  // namespace eigen

namespace plugin {

/// Define the plugin types
/// @{
template <concepts::value V>
struct eigen {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using size_type = algebra::eigen::size_type;
  using transform3D = algebra::eigen::transform3<value_type>;
  using point2D = algebra::eigen::point2<value_type>;
  using point3D = algebra::eigen::point3<value_type>;
  using vector3D = algebra::eigen::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::eigen::matrix_type<value_type, ROWS, COLS>;
};
/// @}

}  // namespace plugin

}  // namespace algebra
