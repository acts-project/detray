/** Detray plugins library, part of the ACTS project
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/generic/impl/generic_matrix.hpp"
#include "detray/algebra/generic/impl/generic_vector.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"

namespace detray::algebra::generic::math {

/// @TODO: Move to algebra plugins
template <concepts::algebra algebra_t>
struct matrix_helper {

    /// Matrix index type
    using index_type = dindex_type<algebra_t>;
    using scalar_type = dscalar<algebra_t>;
    using vector3 = dvector3D<algebra_t>;

    /// Matrix type
    template <index_type ROWS, index_type COLS>
    using matrix_type = dmatrix<algebra_t, ROWS, COLS>;

    /// Cholesky decomposition
    template <concepts::square_matrix matrix_t>
    DETRAY_HOST_DEVICE inline matrix_t cholesky_decomposition(
        const matrix_t& mat) const {

        using index_t = detray::traits::index_t<matrix_t>;
        constexpr index_t N{detray::traits::rank<matrix_t>};

        matrix_t L = detray::algebra::generic::math::zero<matrix_t>();

        // Cholesky–Banachiewicz algorithm
        for (std::size_t i = 0u; i < N; i++) {
            for (std::size_t j = 0u; j <= i; j++) {
                scalar_type sum = 0.f;
                for (std::size_t k = 0u; k < j; k++)
                    sum += getter::element(L, i, k) * getter::element(L, j, k);

                if (i == j) {
                    getter::element(L, i, j) =
                        static_cast<scalar_type>(::detray::algebra::math::sqrt(
                            getter::element(mat, i, i) - sum));
                } else {
                    getter::element(L, i, j) =
                        (1.f / getter::element(L, j, j) *
                         (getter::element(mat, i, j) - sum));
                }
            }
        }

        return L;
    }
};

}  // namespace detray::algebra::generic::math
