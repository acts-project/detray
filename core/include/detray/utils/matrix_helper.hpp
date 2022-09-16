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

template <typename matrix_actor_t>
struct matrix_helper {

    /// Matrix actor
    using matrix_actor = matrix_actor_t;
    /// Size type
    using size_type = typename matrix_actor_t::size_ty;
    /// Scalar type
    using scalar_type = typename matrix_actor_t::scalar_type;
    /// 2D Matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type = typename matrix_actor::template matrix_type<ROWS, COLS>;
    /// Array type
    template <size_type N>
    using array_type = typename matrix_actor::template array_type<N>;
    /// 3-element "vector" type
    using vector3 = array_type<3>;

    /// Column-wise cross product between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> column_wise_cross(const matrix_type<3, 3>& m,
                                               const vector3& v) const {
        matrix_type<3, 3> ret;

        auto m_col0 = matrix_actor().template block<3, 1>(m, 0, 0);
        auto m_col1 = matrix_actor().template block<3, 1>(m, 0, 1);
        auto m_col2 = matrix_actor().template block<3, 1>(m, 0, 2);

        matrix_actor().set_block(ret, vector::cross(m_col0, v), 0, 0);
        matrix_actor().set_block(ret, vector::cross(m_col1, v), 0, 1);
        matrix_actor().set_block(ret, vector::cross(m_col2, v), 0, 2);

        return ret;
    }

    /// Column-wise multiplication between matrix (m) and vector (v)
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> column_wise_multiply(const matrix_type<3, 3>& m,
                                                  const vector3& v) const {
        matrix_type<3, 3> ret;

        for (size_type i = 0; i < 3; i++) {
            for (size_type j = 0; j < 3; j++) {
                matrix_actor().element(ret, j, i) =
                    matrix_actor().element(m, j, i) * v[j];
            }
        }

        return ret;
    }

    /// Cross product matrix
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> cross_matrix(const vector3& v) const {
        matrix_type<3, 3> ret;
        matrix_actor().element(ret, 0, 0) = 0;
        matrix_actor().element(ret, 0, 1) = -v[2];
        matrix_actor().element(ret, 0, 2) = v[1];
        matrix_actor().element(ret, 1, 0) = v[2];
        matrix_actor().element(ret, 1, 1) = 0;
        matrix_actor().element(ret, 1, 2) = -v[0];
        matrix_actor().element(ret, 2, 0) = -v[1];
        matrix_actor().element(ret, 2, 1) = v[0];
        matrix_actor().element(ret, 2, 2) = 0;

        return ret;
    }

    /// Outer product operation
    DETRAY_HOST_DEVICE
    inline matrix_type<3, 3> outer_product(const vector3& v1,
                                           const vector3& v2) const {
        matrix_type<3, 1> m1;
        matrix_actor().element(m1, 0, 0) = v1[0];
        matrix_actor().element(m1, 1, 0) = v1[1];
        matrix_actor().element(m1, 2, 0) = v1[2];

        matrix_type<1, 3> m2;
        matrix_actor().element(m2, 0, 0) = v2[0];
        matrix_actor().element(m2, 0, 1) = v2[1];
        matrix_actor().element(m2, 0, 2) = v2[2];

        return m1 * m2;
    }
};

}  // namespace detray