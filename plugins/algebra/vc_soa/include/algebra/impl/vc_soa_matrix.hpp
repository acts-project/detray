/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/impl/vc_soa_vector.hpp"
#include "detray/algebra/common/concepts.hpp"
#include "detray/algebra/common/matrix.hpp"
#include "detray/algebra/common/qualifiers.hpp"

namespace detray::algebra::vc_soa::math {

using algebra::storage::identity;
using algebra::storage::set_identity;
using algebra::storage::set_zero;
using algebra::storage::transpose;
using algebra::storage::zero;

template <std::size_t ROW, std::size_t COL, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr scalar_t determinant(
    const algebra::storage::matrix<array_t, scalar_t, ROW, COL> &) noexcept {
    // @TODO: Implement
    return scalar_t(0);
}

template <std::size_t ROW, std::size_t COL, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr algebra::storage::matrix<array_t, scalar_t, ROW,
                                                       COL>
inverse(
    const algebra::storage::matrix<array_t, scalar_t, ROW, COL> &m) noexcept {
    // @TODO: Implement
    return m;
}

}  // namespace detray::algebra::vc_soa::math
