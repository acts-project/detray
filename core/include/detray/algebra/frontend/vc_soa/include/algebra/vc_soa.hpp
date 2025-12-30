/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/impl/generic_matrix.hpp"
#include "algebra/math/impl/vc_aos_transform3.hpp"
#include "algebra/math/vc_soa.hpp"
#include "algebra/storage/vc_soa.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

/// @name Operators on @c algebra::storage::vector types
/// @{

using algebra::storage::operator*;
using algebra::storage::operator/;
using algebra::storage::operator-;
using algebra::storage::operator+;

/// @}

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::vc_soa types
/// @{

using vc_soa::storage::block;
using vc_soa::storage::element;
using vc_soa::storage::set_block;
using vc_soa::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc_soa types
/// @{

using vc_soa::math::cross;
using vc_soa::math::dot;
using vc_soa::math::eta;
using vc_soa::math::norm;
using vc_soa::math::normalize;
using vc_soa::math::perp;
using vc_soa::math::phi;
using vc_soa::math::theta;

/// @}

}  // namespace vector

// Produces clash with matrix typedefs in other plugins
namespace matrix {

using vc_soa::math::determinant;
using vc_soa::math::identity;
using vc_soa::math::inverse;
using vc_soa::math::set_identity;
using vc_soa::math::set_zero;
using vc_soa::math::transpose;
using vc_soa::math::zero;

using generic::math::set_inplace_product_left;
using generic::math::set_inplace_product_left_transpose;
using generic::math::set_inplace_product_right;
using generic::math::set_inplace_product_right_transpose;
using generic::math::set_product;
using generic::math::set_product_left_transpose;
using generic::math::set_product_right_transpose;
using generic::math::transposed_product;

}  // namespace matrix

namespace vc_soa {

/// @name Vc based transforms on @c algebra::vc_soa types
/// @{

template <concepts::value T>
using transform3 =
    algebra::vc_aos::math::transform3<algebra::vc_soa::storage_type,
                                      Vc::Vector<T>>;

/// @}

}  // namespace vc_soa

namespace plugin {

/// Define the plugin types
/// @{
template <concepts::value V>
struct vc_soa {
  /// Define scalar precision
  using value_type = V;

  template <concepts::element T>
  using simd = Vc::Vector<T>;

  using boolean = Vc::Mask<V>;

  /// Linear Algebra type definitions
  /// @{
  using scalar = simd<value_type>;
  using size_type = algebra::vc_soa::size_type;
  using transform3D = algebra::vc_soa::transform3<value_type>;
  using point2D = algebra::vc_soa::point2<value_type>;
  using point3D = algebra::vc_soa::point3<value_type>;
  using vector3D = algebra::vc_soa::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::vc_soa::matrix_type<value_type, ROWS, COLS>;
  /// @}
};
/// @}

}  // namespace plugin

}  // namespace algebra
