/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/matrix_helper.hpp"

namespace detray {

/// @brief Helper struct to rotate a vector around a given axis and angle
/// counterclockwisely
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
    /// @returns Get the counterclockwisely-rotated vector
    template <typename vector3_t>
    DETRAY_HOST_DEVICE vector3_t operator()(const vector3_t& v) const {
        return R * v;
    }

    private:
    /// Rotation matrix
    matrix_type<3, 3> R = matrix_operator().template identity<3, 3>();
};

/// @brief Helper struct to perform an euler rotation for a given vector
/// All rotation operations are counterclockwise
template <typename transform3_t>
struct euler_rotation {

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

    // Following the z-x-z convention
    vector3 x{1.0f, 0.f, 0.f};
    vector3 z{0.f, 0.f, 1.f};
    scalar_type alpha{0.f};
    scalar_type beta{0.f};
    scalar_type gamma{0.f};

    /// @returns Get the new x and z axis
    DETRAY_HOST_DEVICE std::pair<vector3, vector3> operator()() const {
        // alpha around z axis
        axis_rotation<transform3_t> axis_rot_alpha(z, alpha);
        auto new_x = axis_rot_alpha(x);

        // beta around x' axis
        axis_rotation<transform3_t> axis_rot_beta(new_x, beta);
        const auto new_z = axis_rot_beta(z);

        axis_rotation<transform3_t> axis_rot_gamma(new_z, gamma);
        new_x = axis_rot_gamma(new_x);

        return std::make_pair(new_x, new_z);
    }

    /// @param v vector to be rotated
    /// @returns Get the rotated vector
    DETRAY_HOST_DEVICE vector3 operator()(const vector3& v0) const {
        // alpha around z axis
        axis_rotation<transform3_t> axis_rot_alpha(z, alpha);
        const auto new_x = axis_rot_alpha(x);
        const auto v1 = axis_rot_alpha(v0);

        // beta around x' axis
        axis_rotation<transform3_t> axis_rot_beta(new_x, beta);
        const auto new_z = axis_rot_beta(z);
        const auto v2 = axis_rot_beta(v1);

        // gamma around z' axis
        axis_rotation<transform3_t> axis_rot_gamma(new_z, gamma);
        const auto v3 = axis_rot_gamma(v2);

        return v3;
    }
};

}  // namespace detray
