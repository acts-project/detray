/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/array_cmath.hpp"

#define __plugin algebra::array
#define ALGEBRA_PLUGIN array

namespace detray {

// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
using transform3D = algebra::array::transform3<scalar>;
using point3D = algebra::array::point3<scalar>;
using vector3D = algebra::array::vector3<scalar>;
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
