#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <math.h>
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

/// @brief Calculates the euler angles of a rotation matrix.
///
/// @param matrix The rotation matrix.
///
/// @returns an array containing the euler angles around the x, y, and z axis (in that order).
template <typename matrix_t>
inline std::array<detray::scalar, 3> rotation_matrix_to_euler_angles(
    const matrix_t& matrix) {
    float a = std::sqrt(matrix[0][0] * matrix[0][0] + matrix[1][0] * matrix[1][0]);
    // Checking if it is singular.
    if (a < 1e-6)
        return {std::atan2(-matrix[1][2], matrix[1][1]), std::atan2(-matrix[2][0], a),
                0};

    return {std::atan2(matrix[2][1], matrix[2][2]), std::atan2(-matrix[2][0], a),
            std::atan2(matrix[1][0], matrix[0][0])};
}

/// @brief Calculates the detray point3 as an actsvg point.
///
/// @param d_point The detray point3.
///
/// @returns an actsvg point3.
template <std::size_t dim, typename point_t>
inline std::array<actsvg::scalar, dim> convert_point(
    const point_t& d_point) {
    std::array<actsvg::scalar, dim> ret;
    for (std::size_t i = 0; i < ret.size(); i++){
        ret[i] = static_cast<actsvg::scalar>(d_point[i]);
    }
    return ret;
}

/// @brief Calculates the equivalent actsvg transform of a detray transform.
///
/// @param d_transform The detray transform.
///
/// @returns an actsvg transform.
template <typename transform_t>
inline auto convert_transform(const transform_t& d_transform) {
    auto translation = d_transform.translation();
    auto euler_angles =
        rotation_matrix_to_euler_angles<>(d_transform.rotation());

    // TODO: skew and scale

    auto ret = actsvg::style::transform();
    constexpr auto rad_to_deg = 180.0 / 3.14;
    ret._tr = {static_cast<actsvg::scalar>(translation[0]),
               static_cast<actsvg::scalar>(translation[1])};
    ret._rot = {static_cast<actsvg::scalar>(euler_angles[0] * rad_to_deg),
                static_cast<actsvg::scalar>(euler_angles[1] * rad_to_deg),
                static_cast<actsvg::scalar>(euler_angles[2] * rad_to_deg)};

    return ret;
}
}  // namespace detray::actsvg_visualization