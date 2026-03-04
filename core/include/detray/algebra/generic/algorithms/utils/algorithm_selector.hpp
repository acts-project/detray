/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/generic/algorithms/matrix/determinant/hard_coded.hpp"
#include "detray/algebra/generic/algorithms/matrix/determinant/partial_pivot_lud.hpp"
#include "detray/algebra/generic/algorithms/matrix/inverse/hard_coded.hpp"
#include "detray/algebra/generic/algorithms/matrix/inverse/partial_pivot_lud.hpp"

namespace detray::algebra::generic {

/// Get the type of determinant algorithm according to matrix dimension
/// @{

// Default algorithm
template <std::size_t N, typename... Args>
struct determinant_selector {
    using type = matrix::determinant::partial_pivot_lud<Args...>;
};

/// Always use hard coded implementation for very small matrices
template <typename... Args>
struct determinant_selector<2, Args...> {
    using type = matrix::determinant::hard_coded<Args...>;
};

/// @tparam M matrix type
template <concepts::square_matrix M>
using determinant_t =
    typename determinant_selector<detray::traits::rank<M>, M>::type;
/// @}

/// Get the type of inversion algorithm according to matrix dimension
/// @{
template <std::size_t N, typename... Args>
struct inversion_selector {
    using type = matrix::inverse::partial_pivot_lud<Args...>;
};

/// Always use hard coded implementation for very small matrices
template <typename... Args>
struct inversion_selector<2, Args...> {
    using type = matrix::inverse::hard_coded<Args...>;
};

/// @tparam M matrix type
template <concepts::square_matrix M>
using inversion_t =
    typename inversion_selector<detray::traits::rank<M>, M>::type;
/// @}

}  // namespace detray::algebra::generic
