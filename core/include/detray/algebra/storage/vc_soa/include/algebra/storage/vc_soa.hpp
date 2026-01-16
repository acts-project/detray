/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/storage/impl/vc_aos_approximately_equal.hpp"
#include "algebra/storage/impl/vc_soa_casts.hpp"
#include "algebra/storage/impl/vc_soa_getter.hpp"
#include "algebra/storage/matrix.hpp"
#include "algebra/storage/vector.hpp"
#include "algebra/type_traits.hpp"

// System include(s).
#include <array>
#include <cstddef>

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace algebra {

namespace vc_soa {

/// Size type for Vc storage model
using size_type = std::size_t;
/// Value type in a linear algebra vector: SoA layout
template <concepts::value T>
using value_type = T;
/// Scalar type in a linear algebra vector: SoA layout
template <concepts::value T>
using scalar_type = Vc::Vector<T>;
/// Array type used to store Vc::Vectors or matrix columns
template <concepts::simd_scalar T, size_type N>
using storage_type = std::array<T, N>;
/// Vector type used in the Vc SoA storage model
template <concepts::value T, std::size_t N>
using vector_type = algebra::storage::vector<N, Vc::Vector<T>, storage_type>;
/// Matrix type used in the Vc SoA storage model
template <concepts::value T, size_type ROWS, size_type COLS>
using matrix_type =
    algebra::storage::matrix<storage_type, Vc::Vector<T>, ROWS, COLS>;

/// 2-element "vector" type, using @c Vc::Vector in every element
template <concepts::value T>
using vector2 = vector_type<T, 2>;
/// Point in 2D space, using @c Vc::Vector in every element
template <concepts::value T>
using point2 = vector2<T>;
/// 3-element "vector" type, using @c Vc::Vector in every element
template <concepts::value T>
using vector3 = vector_type<T, 3>;
/// Point in 3D space, using @c Vc::Vector in every element
template <concepts::value T>
using point3 = vector3<T>;
/// 6-element "vector" type, using @c Vc::Vector in every element
template <concepts::value T>
using vector6 = vector_type<T, 6>;
/// 8-element "vector" type, using @c Vc::Vector in every element
template <concepts::value T>
using vector8 = vector_type<T, 8>;

/// Element Getter
using element_getter = algebra::storage::element_getter;
/// Block Getter
using block_getter = algebra::storage::block_getter;

}  // namespace vc_soa

ALGEBRA_PLUGINS_DEFINE_TYPE_TRAITS(vc_soa)

namespace traits {

/// Make sure the simd scalar type is recognized correctly
/// @{
template <concepts::value T, std::size_t ROWS, std::size_t COLS>
struct scalar<
    algebra::storage::matrix<vc_soa::storage_type, Vc::Vector<T>, ROWS, COLS>> {
  using type = Vc::Vector<T>;
};

template <concepts::value T, std::size_t N>
struct scalar<
    algebra::storage::vector<N, Vc::Vector<T>, vc_soa::storage_type>> {
  using type = Vc::Vector<T>;
};
/// @}

/// Get the single value type from the simd scalar type
/// @{
template <concepts::value T>
struct value<Vc::Vector<T>> {
  using type = T;
};
/// @}

// Vector and storage types are different
template <concepts::simd_scalar T, auto N>
struct dimensions<vc_soa::storage_type<T, N>> {

  using size_type = vc_soa::size_type;

  static constexpr size_type dim{1};
  static constexpr size_type rows{N};
  static constexpr size_type columns{1};
};

}  // namespace traits

}  // namespace algebra
