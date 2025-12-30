/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/concepts.hpp"
#include "algebra/storage/impl/smatrix_getter.hpp"
#include "algebra/type_traits.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>
#include <Math/SVector.h>

// System include(s).
#include <cstddef>

namespace algebra {

namespace smatrix {

/// size type for SMatrix storage model
using size_type = unsigned int;
/// Value type for SMatrix storage model
template <concepts::value T>
using value_type = T;
/// Scalar type for SMatrix storage model
template <concepts::value T>
using scalar_type = T;
/// Array type used in the SMatrix storage model
template <concepts::scalar T, size_type N>
using storage_type = ROOT::Math::SVector<T, N>;
/// Vector type used in the SMatrix storage model
template <concepts::scalar T, size_type N>
using vector_type = storage_type<T, N>;
/// Matrix type used in the SMatrix storage model
template <concepts::scalar T, size_type ROWS, size_type COLS>
using matrix_type = ROOT::Math::SMatrix<T, ROWS, COLS>;

/// 3-element "vector" type, using @c ROOT::Math::SVector
template <concepts::scalar T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c ROOT::Math::SVector
template <concepts::scalar T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c ROOT::Math::SVector
template <concepts::scalar T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c ROOT::Math::SVector
template <concepts::scalar T>
using point2 = vector2<T>;

/// Element Getter
using element_getter = smatrix::storage::element_getter;
/// Block Getter
using block_getter = smatrix::storage::block_getter;

}  // namespace smatrix

ALGEBRA_PLUGINS_DEFINE_TYPE_TRAITS(smatrix)

}  // namespace algebra
