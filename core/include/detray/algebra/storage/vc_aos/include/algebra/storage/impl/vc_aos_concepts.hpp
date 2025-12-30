/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/vector.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s).
#include <concepts>

namespace algebra::concepts {

/// Vc SIMD types
template <typename T>
concept vc_simd_vector = Vc::is_simd_vector<T>::value;

/// Vc SIMD types
template <typename T>
concept simd_storage_vector = algebra::detail::is_storage_vector_v<T>;

/// Vc AoS vector
template <typename T>
concept vc_aos_vector = (vc_simd_vector<T> || simd_storage_vector<T>);

}  // namespace algebra::concepts
