/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/definitions/qualifiers.hpp"

namespace detray {

template <typename scalar_t>
struct vector_helpers {

    using vector3 = __plugin::vector3<scalar_t>;

    using size_type = __plugin::size_type;

    template <size_type ROWS, size_type COLS>
    using matrix_type = __plugin::matrix_type<scalar_t, ROWS, COLS>;

    using matrix_operator = standard_matrix_operator<scalar_t>;

    /// Column-wise cross product between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
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

    /// Column-wise dot product between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> dot(const matrix_type<3, 3>& m,
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