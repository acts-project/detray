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

// Define namespace(s)
namespace getter = algebra::getter;
namespace vector = algebra::vector;
namespace matrix = algebra::matrix;

// Define matrix/vector operator
template <typename scalar_t>
using standard_matrix_operator = matrix::actor<scalar_t>;

}  // namespace detray
