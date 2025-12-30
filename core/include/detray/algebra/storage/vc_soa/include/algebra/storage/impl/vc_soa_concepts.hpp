/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/impl/vc_aos_concepts.hpp"

// System include(s).
#include <concepts>

namespace algebra::concepts {

/// Vc SoA vector
template <typename T>
concept vc_soa_vector = (simd_storage_vector<T> ||
                         vc_simd_vector<typename T::value_type>);

}  // namespace algebra::concepts
