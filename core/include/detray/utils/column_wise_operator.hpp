/** Detray plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"

namespace detray {

template <typename matrix_operator_t>
struct column_wise_operator {

    /// Matrix actor
    using matrix_operator = matrix_operator_t;
    /// Size type
    using size_type = typename matrix_operator_t::size_ty;
    /// Scalar type
    using scalar_type = typename matrix_operator_t::scalar_type;
    /// 2D Matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    /// Array type
    template <size_type N>
    using array_type = typename matrix_operator::template array_type<N>;
    /// 3-element "vector" type
    using vector3 = array_type<3>;

    /// Column-wise cross product between matrix (m) and vector (v)
    ALGEBRA_HOST_DEVICE
    inline matrix_type<3, 3> cross(const matrix_type<3, 3>& m,
                                   const vector3& v) const {
        matrix_type<3, 3> ret;

        auto m_col0 = matrix_operator().template block<3, 1>(m, 0, 0);
        auto m_col1 = matrix_operator().template block<3, 1>(m, 0, 1);
        auto m_col2 = matrix_operator().template block<3, 1>(m, 0, 2);

        matrix_operator().set_block(ret, vector::cross(m_col0, v), 0, 0);
        matrix_operator().set_block(ret, vector::cross(m_col1, v), 0, 1);
        matrix_operator().set_block(ret, vector::cross(m_col2, v), 0, 2);

        return ret;
    }

    /// Column-wise multiplication between matrix (m) and vector (v)
    ALGEBRA_HOST_DEVICE
    inline matrix_type<3, 3> multiply(const matrix_type<3, 3>& m,
                                      const vector3& v) const {
        matrix_type<3, 3> ret;

        for (size_type i = 0; i < 3; i++) {
            for (size_type j = 0; j < 3; j++) {
                matrix_operator().element(ret, j, i) =
                    matrix_operator().element(m, j, i) * v[j];
            }
        }

        return ret;
    }
};

}  // namespace detray