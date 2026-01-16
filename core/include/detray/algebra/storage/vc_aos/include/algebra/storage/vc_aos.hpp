/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/storage/impl/vc_aos_approximately_equal.hpp"
#include "algebra/storage/impl/vc_aos_concepts.hpp"
#include "algebra/storage/impl/vc_aos_getter.hpp"
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

namespace vc_aos {

/// Size type for Vc storage model
using size_type = std::size_t;
/// Value type for Vc storage model
template <concepts::value T>
using value_type = T;
/// Scalar type for Vc storage model
template <concepts::value T>
using scalar_type = T;
/// Array type used to store Vc::Vectors
template <concepts::value T, size_type N>
using storage_type = Vc::SimdArray<T, N>;
/// Value type in a linear algebra vector: AoS layout
template <concepts::value T>
using value_type = T;
/// Scalar type in a linear algebra vector: AoS layout
template <concepts::value T>
using scalar_type = T;
/// Vector type used in the Vc AoS storage model
template <concepts::value T, std::size_t N>
using vector_type = algebra::storage::vector<N, T, storage_type>;
/// Matrix type used in the Vc AoS storage model
template <concepts::value T, size_type ROWS, size_type COLS>
using matrix_type = algebra::storage::matrix<storage_type, T, ROWS, COLS>;

/// 2-element "vector" type, using @c Vc::SimdArray
template <concepts::value T>
using vector2 = vector_type<T, 2>;
/// Point in 2D space, using @c Vc::SimdArray
template <concepts::value T>
using point2 = vector2<T>;
/// 3-element "vector" type, using @c Vc::SimdArray
template <concepts::value T>
using vector3 = vector_type<T, 3>;
/// Point in 3D space, using @c Vc::SimdArray
template <concepts::value T>
using point3 = vector3<T>;
/// 6-element "vector" type, using @c Vc::SimdArray
template <concepts::value T>
using vector6 = vector_type<T, 6>;
/// 8-element "vector" type, using @c Vc::SimdArray
template <concepts::value T>
using vector8 = vector_type<T, 8>;

/// Element Getter
using element_getter = algebra::storage::element_getter;
/// Block Getter
using block_getter = algebra::storage::block_getter;

}  // namespace vc_aos

ALGEBRA_PLUGINS_DEFINE_TYPE_TRAITS(vc_aos)

namespace traits {

template <typename T, auto N>
struct index<vc_aos::storage_type<T, N>> {
  using size_type = vc_aos::size_type;
};

template <typename T, auto N>
struct index<Vc_1::Vector<T, Vc_1::simd_abi::fixed_size<N>>> {
  using size_type = vc_aos::size_type;
};

template <typename T, auto N>
struct value<vc_aos::storage_type<T, N>> {
  using type = T;
};

template <typename T, auto N>
struct value<Vc_1::Vector<T, Vc_1::simd_abi::fixed_size<N>>> {
  using type = T;
};

// Vector and storage types are different
template <typename T, auto N>
struct dimensions<vc_aos::storage_type<T, N>> {

  using size_type = vc_aos::size_type;

  static constexpr size_type dim{1};
  static constexpr size_type rows{N};
  static constexpr size_type columns{1};
};

// Vector and storage types are different
template <typename T, auto N>
struct dimensions<Vc_1::Vector<T, Vc_1::simd_abi::fixed_size<N>>> {

  using size_type = vc_aos::size_type;

  static constexpr size_type dim{1};
  static constexpr size_type rows{N};
  static constexpr size_type columns{1};
};

}  // namespace traits

}  // namespace algebra
