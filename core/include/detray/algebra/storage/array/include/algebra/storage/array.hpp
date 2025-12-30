/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/concepts.hpp"
#include "algebra/storage/impl/cmath_getter.hpp"
#include "algebra/type_traits.hpp"

// System include(s).
#include <array>
#include <cstddef>

namespace algebra {

namespace array {

/// size type for Array storage model
using size_type = std::size_t;
/// Value type for Array storage model
template <concepts::value T>
using value_type = T;
/// Scalar type for Array storage model
template <concepts::value T>
using scalar_type = T;
/// Array type used in the Array storage model
template <typename T, size_type N>
using storage_type = std::array<T, N>;
/// Vector type used in the Array storage model
template <concepts::scalar T, std::size_t N>
using vector_type = storage_type<T, N>;
/// Matrix type used in the Array storage model
template <concepts::scalar T, size_type ROWS, size_type COLS>
using matrix_type = storage_type<storage_type<T, ROWS>, COLS>;

/// 3-element "vector" type, using @c std::array
template <concepts::scalar T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c std::array
template <concepts::scalar T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c std::array
template <concepts::scalar T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c std::array
template <concepts::scalar T>
using point2 = vector2<T>;

/// Element Getter
using element_getter = cmath::storage::element_getter;
/// Block Getter
using block_getter = cmath::storage::block_getter;

}  // namespace array

ALGEBRA_PLUGINS_DEFINE_TYPE_TRAITS(array)

}  // namespace algebra
