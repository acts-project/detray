/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/matrix_helper.hpp"

namespace detray {

template <typename transform3_t>
struct axis_rotation {
    using point3 = typename transform3_t::point3;
    using vector3 = typename transform3_t::vector3;
    using scalar_type = typename transform3_t::scalar_type;
    using matrix_operator = typename transform3_t::matrix_actor;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using mat_helper = matrix_helper<matrix_operator>;

    DETRAY_HOST_DEVICE
    axis_rotation(const vector3& axis, const scalar_type theta) {
        scalar_type cos_theta = math_ns::cos(theta);

        matrix_type<3, 3> I = matrix_operator().template identity<3, 3>();
        matrix_type<3, 3> axis_cross = mat_helper().cross_matrix(axis);
        matrix_type<3, 3> axis_outer = mat_helper().outer_product(axis, axis);

        R = cos_theta * I + std::sin(theta) * axis_cross +
            (1 - cos_theta) * axis_outer;
    }

    template <typename vector3_t>
    DETRAY_HOST_DEVICE vector3_t operator()(const vector3_t& v) {
        return R * v;
    }

    matrix_type<3, 3> R = matrix_operator().template identity<3, 3>();
};

}  // namespace detray