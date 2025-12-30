/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/storage/impl/fastor_getter.hpp"
#include "algebra/storage/impl/fastor_matrix.hpp"
#include "algebra/type_traits.hpp"

// System include(s).
#include <cstddef>

namespace algebra {

namespace fastor {

/// size type for Fastor storage model
using size_type = std::size_t;
/// Value type for Fastor storage model
template <concepts::value T>
using value_type = T;
/// Scalar type for Fastor storage model
template <concepts::value T>
using scalar_type = T;
/// Array type used in the Fastor storage model
template <concepts::scalar T, size_type N>
using storage_type = Fastor::Tensor<T, N>;
/// Vector type used in the Fastor storage model
template <concepts::scalar T, size_type N>
using vector_type = storage_type<T, N>;
/// Matrix type used in the Fastor storage model
template <concepts::scalar T, size_type ROWS, size_type COLS>
using matrix_type = algebra::fastor::Matrix<T, ROWS, COLS>;

/// 3-element "vector" type, using @c Fastor::Tensor
template <concepts::scalar T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c Fastor::Tensor
template <concepts::scalar T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c Fastor::Tensor
template <concepts::scalar T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c Fastor::Tensor
template <concepts::scalar T>
using point2 = vector2<T>;

/// Element Getter
using element_getter = fastor::storage::element_getter;
/// Block Getter
using block_getter = fastor::storage::block_getter;

}  // namespace fastor

ALGEBRA_PLUGINS_DEFINE_TYPE_TRAITS(fastor)

}  // namespace algebra
