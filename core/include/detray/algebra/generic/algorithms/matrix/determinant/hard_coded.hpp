/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/qualifiers.hpp"
#include "detray/algebra/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::algebra::generic::matrix::determinant {

/// "Determinant getter", assuming a N X N matrix
template <concepts::square_matrix matrix_t>
struct hard_coded {

    using scalar_t = algebra::traits::value_t<matrix_t>;
    using index_t = algebra::traits::index_t<matrix_t>;

    /// Function (object) used for accessing a matrix element
    using element_getter_t = algebra::traits::element_getter_t<matrix_t>;

    // 2 X 2 matrix determinant
    template <typename M = matrix_t>
        requires(algebra::traits::rank<M> == 2)
    ALGEBRA_HOST_DEVICE constexpr scalar_t operator()(const matrix_t &m) const {

        constexpr element_getter_t elem{};

        return elem(m, 0, 0) * elem(m, 1, 1) - elem(m, 0, 1) * elem(m, 1, 0);
    }

    // 4 X 4 matrix determinant
    template <typename M = matrix_t>
        requires(algebra::traits::rank<M> == 4)
    ALGEBRA_HOST_DEVICE constexpr scalar_t operator()(const matrix_t &m) const {

        constexpr element_getter_t elem{};

        return elem(m, 0, 3) * elem(m, 1, 2) * elem(m, 2, 1) * elem(m, 3, 0) -
               elem(m, 0, 2) * elem(m, 1, 3) * elem(m, 2, 1) * elem(m, 3, 0) -
               elem(m, 0, 3) * elem(m, 1, 1) * elem(m, 2, 2) * elem(m, 3, 0) +
               elem(m, 0, 1) * elem(m, 1, 3) * elem(m, 2, 2) * elem(m, 3, 0) +
               elem(m, 0, 2) * elem(m, 1, 1) * elem(m, 2, 3) * elem(m, 3, 0) -
               elem(m, 0, 1) * elem(m, 1, 2) * elem(m, 2, 3) * elem(m, 3, 0) -
               elem(m, 0, 3) * elem(m, 1, 2) * elem(m, 2, 0) * elem(m, 3, 1) +
               elem(m, 0, 2) * elem(m, 1, 3) * elem(m, 2, 0) * elem(m, 3, 1) +
               elem(m, 0, 3) * elem(m, 1, 0) * elem(m, 2, 2) * elem(m, 3, 1) -
               elem(m, 0, 0) * elem(m, 1, 3) * elem(m, 2, 2) * elem(m, 3, 1) -
               elem(m, 0, 2) * elem(m, 1, 0) * elem(m, 2, 3) * elem(m, 3, 1) +
               elem(m, 0, 0) * elem(m, 1, 2) * elem(m, 2, 3) * elem(m, 3, 1) +
               elem(m, 0, 3) * elem(m, 1, 1) * elem(m, 2, 0) * elem(m, 3, 2) -
               elem(m, 0, 1) * elem(m, 1, 3) * elem(m, 2, 0) * elem(m, 3, 2) -
               elem(m, 0, 3) * elem(m, 1, 0) * elem(m, 2, 1) * elem(m, 3, 2) +
               elem(m, 0, 0) * elem(m, 1, 3) * elem(m, 2, 1) * elem(m, 3, 2) +
               elem(m, 0, 1) * elem(m, 1, 0) * elem(m, 2, 3) * elem(m, 3, 2) -
               elem(m, 0, 0) * elem(m, 1, 1) * elem(m, 2, 3) * elem(m, 3, 2) -
               elem(m, 0, 2) * elem(m, 1, 1) * elem(m, 2, 0) * elem(m, 3, 3) +
               elem(m, 0, 1) * elem(m, 1, 2) * elem(m, 2, 0) * elem(m, 3, 3) +
               elem(m, 0, 2) * elem(m, 1, 0) * elem(m, 2, 1) * elem(m, 3, 3) -
               elem(m, 0, 0) * elem(m, 1, 2) * elem(m, 2, 1) * elem(m, 3, 3) -
               elem(m, 0, 1) * elem(m, 1, 0) * elem(m, 2, 2) * elem(m, 3, 3) +
               elem(m, 0, 0) * elem(m, 1, 1) * elem(m, 2, 2) * elem(m, 3, 3);
    }
};

}  // namespace detray::algebra::generic::matrix::determinant
