/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/vc_aos.hpp"

namespace detray {

/// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct vc_aos {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using transform3D = algebra::vc_aos::transform3<value_type>;
    using point2D = algebra::vc_aos::point2<value_type>;
    using point3D = algebra::vc_aos::point3<value_type>;
    using vector3D = algebra::vc_aos::vector3<value_type>;
};
/// @}

namespace vector {

using algebra::vc_aos::math::cross;
using algebra::vc_aos::math::dot;
using algebra::vc_aos::math::normalize;

}  // namespace vector

namespace getter {

using algebra::vc_aos::math::eta;
using algebra::vc_aos::math::norm;
using algebra::vc_aos::math::perp;
using algebra::vc_aos::math::phi;
using algebra::vc_aos::math::theta;

using algebra::cmath::element;

/// Function extracting a slice from the matrix used by
/// @c algebra::vc_aos::transform3<float>
template <std::size_t SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline algebra::vc_aos::vector3<float> vector(
    const algebra::vc_aos::transform3<float>::matrix44& m,
    std::size_t
#ifndef NDEBUG
        row
#endif  // not NDEBUG
    ,
    std::size_t col) {

    assert(row == 0);
    assert(col < 4);
    switch (col) {
        case 0:
            return m.x;
        case 1:
            return m.y;
        case 2:
            return m.z;
        case 3:
            return m.t;
        default:
            return m.x;
    }
}

/// Function extracting a slice from the matrix used by
/// @c algebra::vc_aos::transform3<double>
template <std::size_t SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline algebra::vc_aos::vector3<double> vector(
    const algebra::vc_aos::transform3<double>::matrix44& m,
    std::size_t
#ifndef NDEBUG
        row
#endif  // not NDEBUG
    ,
    std::size_t col) {

    assert(row == 0);
    assert(col < 4);
    switch (col) {
        case 0:
            return m.x;
        case 1:
            return m.y;
        case 2:
            return m.z;
        case 3:
            return m.t;
        default:
            return m.x;
    }
}

}  // namespace getter

}  // namespace detray
