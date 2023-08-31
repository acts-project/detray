/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/vc_cmath.hpp"

#define __plugin algebra::vc
#define ALGEBRA_PLUGIN vc_array

namespace detray {

/// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
using transform3D = algebra::vc::transform3<scalar>;
using point3D = algebra::vc::point3<scalar>;
using vector3D = algebra::vc::vector3<scalar>;
/// @}

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct vc {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using transform3D = algebra::vc::transform3<value_type>;
    using point2D = algebra::vc::point2<value_type>;
    using point3D = algebra::vc::point3<value_type>;
    using vector3D = algebra::vc::vector3<value_type>;
};
/// @}

// Define namespace(s)
namespace getter = algebra::getter;
namespace vector = algebra::vector;
namespace matrix = algebra::matrix;

// Define matrix/vector operator
template <typename scalar_t>
using standard_matrix_operator =
    matrix::actor<scalar_t, matrix::determinant::preset0<scalar_t>,
                  matrix::inverse::preset0<scalar_t>>;

}  // namespace detray
