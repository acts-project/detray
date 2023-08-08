// Not used currently in visualization implementation i.e. detray-actsvg
// conversion.

#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"

// System include(s)
#include <math.h>

#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization::proto::utils {

/// @brief Calculates the euler angles of a rotation matrix.
///
/// @param matrix The rotation matrix.
///
/// @returns An array containing the euler angles around the x, y, and z axis
/// (in that order).
template <typename matrix_t>
inline std::array<detray::scalar, 3> mat_to_euler(const matrix_t& matrix) {
    float a =
        std::sqrt(matrix[0][0] * matrix[0][0] + matrix[1][0] * matrix[1][0]);
    // Checking if it is singular.
    if (a < 1e-6) {
        return {std::atan2(-matrix[1][2], matrix[1][1]),
                std::atan2(-matrix[2][0], a), 0};
    }
    return {std::atan2(matrix[2][1], matrix[2][2]),
            std::atan2(-matrix[2][0], a),
            std::atan2(matrix[1][0], matrix[0][0])};
}

/// @brief Calculates the detray point3 as an actsvg point.
///
/// @param d_point The detray point3.
///
/// @returns An actsvg point3.
template <std::size_t dim, typename point_t>
inline std::array<actsvg::scalar, dim> convert_point(const point_t& d_point) {
    std::array<actsvg::scalar, dim> ret;
    for (std::size_t i = 0; i < ret.size(); i++) {
        ret[i] = static_cast<actsvg::scalar>(d_point[i]);
    }
    return ret;
}

/// @brief Calculates the equivalent actsvg transform of a detray transform.
///
/// @param d_transform The detray transform.
///
/// @note Does not apply skew, scale or rotation around the x or y axis.
///
/// @returns An actsvg transform.
template <typename transform_t>
inline auto to_actsvg_transform(const transform_t& d_transform) {
    auto translation = d_transform.translation();
    auto euler_angles = rotation_matrix_to_euler_angles(d_transform.rotation());

    auto ret = actsvg::style::transform();

    // The translate(<x> [<y>]) transform function moves the object by x and y.
    ret._tr = {static_cast<actsvg::scalar>(translation[0]),
               static_cast<actsvg::scalar>(translation[1])};

    // The rotate(<a> [<x> <y>]) transform function specifies a rotation by a
    // degrees about a given point which is (0,0) here.
    ret._rot = {static_cast<actsvg::scalar>(
                    euler_angles[2] * detray::unit<detray::scalar>::rad_to_deg),
                static_cast<actsvg::scalar>(0), static_cast<actsvg::scalar>(0)};

    return ret;
}
}  // namespace detray::actsvg_visualization::proto::utils