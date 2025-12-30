/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace algebra::generic::matrix::determinant {

/// "Determinant getter", assuming a N X N matrix
template <concepts::square_matrix matrix_t, class element_getter_t>
struct cofactor {

  using scalar_type = algebra::traits::value_t<matrix_t>;
  using size_type = algebra::traits::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  ALGEBRA_HOST_DEVICE constexpr scalar_type operator()(
      const matrix_t &m) const {
    return determinant_getter_helper<algebra::traits::rank<matrix_t>>()(m);
  }

  template <size_type N, typename Enable = void>
  struct determinant_getter_helper;

  template <size_type N>
  struct determinant_getter_helper<N, typename std::enable_if_t<N == 1>> {
    template <class input_matrix_type>
    ALGEBRA_HOST_DEVICE constexpr scalar_type operator()(
        const input_matrix_type &m) const {
      return element_getter()(m, 0, 0);
    }
  };

  template <size_type N>
  struct determinant_getter_helper<N, typename std::enable_if_t<N != 1>> {

    template <class input_matrix_type>
    ALGEBRA_HOST_DEVICE constexpr scalar_type operator()(
        const input_matrix_type &m) const {

      scalar_type D = 0;

      // To store cofactors
      matrix_t temp;

      // To store sign multiplier
      int sign = 1;

      // Iterate for each element of first row
      for (size_type col = 0; col < N; col++) {
        // Getting Cofactor of A[0][f]
        this->get_cofactor(m, temp, size_type(0), col);
        D += sign * element_getter()(m, 0, col) *
             determinant_getter_helper<N - 1>()(temp);

        // terms are to be added with alternate sign
        sign = -sign;
      }

      return D;
    }

    template <class input_matrix_type>
    ALGEBRA_HOST_DEVICE constexpr void get_cofactor(const input_matrix_type &m,
                                                    matrix_t &temp, size_type p,
                                                    size_type q) const {

      size_type i = 0;
      size_type j = 0;

      // Looping for each element of the matrix
      for (size_type row = 0; row < N; row++) {
        for (size_type col = 0; col < N; col++) {
          //  Copying into temporary matrix only those element
          //  which are not in given row and column
          if (row != p && col != q) {
            element_getter()(temp, i, j++) = element_getter()(m, row, col);

            // Row is filled, so increase row index and
            // reset col index
            if (j == N - 1) {
              j = 0;
              i++;
            }
          }
        }
      }
    }
  };
};

}  // namespace algebra::generic::matrix::determinant
