/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/matrix_helper.hpp"

namespace detray {

/// @brief Helper struct to rotate a vector around a given axis and angle
template <typename transform3_t>
struct axis_rotation {

    public:
    using point3 = typename transform3_t::point3;
    using vector3 = typename transform3_t::vector3;
    using scalar_type = typename transform3_t::scalar_type;
    using matrix_operator = typename transform3_t::matrix_actor;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using mat_helper = matrix_helper<matrix_operator>;

    /// @brief Constructor for axis rotation
    ///
    /// @param axis rotation axis
    /// @param theta rotation angle
    DETRAY_HOST_DEVICE
    axis_rotation(const vector3& axis, const scalar_type theta) {
        // normalize the axis
        const auto U = vector::normalize(axis);

        scalar_type cos_theta{math::cos(theta)};

        matrix_type<3, 3> I = matrix_operator().template identity<3, 3>();
        matrix_type<3, 3> axis_cross = mat_helper().cross_matrix(U);
        matrix_type<3, 3> axis_outer = mat_helper().outer_product(U, U);

        R = cos_theta * I + math::sin(theta) * axis_cross +
            (1.f - cos_theta) * axis_outer;
    }

    /// @param v vector to be rotated
    /// @returns Get the rotated vector
    template <typename vector3_t>
    DETRAY_HOST_DEVICE vector3_t operator()(const vector3_t& v) const {
        return R * v;
    }

    private:
    /// Rotation matrix
    matrix_type<3, 3> R = matrix_operator().template identity<3, 3>();
};

}  // namespace detray
