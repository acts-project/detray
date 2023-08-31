/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/smatrix_smatrix.hpp"

#define __plugin algebra::smatrix
#define ALGEBRA_PLUGIN smatrix

namespace detray {

// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
using transform3D = algebra::smatrix::transform3<scalar>;
using point3D = algebra::smatrix::point3<scalar>;
using vector3D = algebra::smatrix::vector3<scalar>;
/// @}

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct smatrix {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using transform3D = algebra::smatrix::transform3<value_type>;
    using point2D = algebra::smatrix::point2<value_type>;
    using point3D = algebra::smatrix::point3<value_type>;
    using vector3D = algebra::smatrix::vector3<value_type>;
};
/// @}

// Define namespace(s)
namespace getter = algebra::getter;
namespace vector = algebra::vector;
namespace matrix = algebra::matrix;

// Define matrix/vector operator
template <typename scalar_t>
using standard_matrix_operator = matrix::actor<scalar_t>;

}  // namespace detray
