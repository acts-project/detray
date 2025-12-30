/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/impl/vc_soa_vector.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/matrix.hpp"

namespace algebra::vc_soa::math {

using storage::identity;
using storage::set_identity;
using storage::set_zero;
using storage::transpose;
using storage::zero;

template <std::size_t ROW, std::size_t COL, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr scalar_t determinant(
    const storage::matrix<array_t, scalar_t, ROW, COL> &) noexcept {
  // @TODO: Implement
  return scalar_t(0);
}

template <std::size_t ROW, std::size_t COL, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr storage::matrix<array_t, scalar_t, ROW, COL>
inverse(const storage::matrix<array_t, scalar_t, ROW, COL> &m) noexcept {
  // @TODO: Implement
  return m;
}

}  // namespace algebra::vc_soa::math
