/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/eigen_eigen.hpp"

#define __plugin algebra::eigen
#define ALGEBRA_PLUGIN eigen

namespace detray {

// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

// Define namespace(s)
namespace getter = algebra::getter;
namespace vector = algebra::vector;
namespace matrix = algebra::matrix;

// Define matrix/vector operator
template <typename scalar_t>
using standard_matrix_operator = matrix::actor<scalar_t>;

}  // namespace detray