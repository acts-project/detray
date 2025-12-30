/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/generic.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/matrix.hpp"

namespace algebra::vc_aos::math {

using storage::identity;
using storage::set_identity;
using storage::set_zero;
using storage::transpose;
using storage::zero;

/// @returns the determinant
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr value_t determinant(
    const storage::matrix<array_t, value_t, N, N> &m) noexcept {
  return algebra::generic::math::determinant(m);
}

/// @returns the inverse
template <std::size_t ROW, std::size_t COL, concepts::value value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr storage::matrix<array_t, value_t, COL, ROW>
inverse(const storage::matrix<array_t, value_t, ROW, COL> &m) noexcept {
  return algebra::generic::math::inverse(m);
}

/// @returns the transpose
template <std::size_t ROW, std::size_t COL, concepts::value value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr storage::matrix<array_t, value_t, COL, ROW>
transpose(const storage::matrix<array_t, value_t, ROW, COL> &m) noexcept {
  return algebra::generic::math::transpose(m);
}

}  // namespace algebra::vc_aos::math
