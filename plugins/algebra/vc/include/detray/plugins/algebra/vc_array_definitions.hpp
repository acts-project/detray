/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra-Plugins include
#include "algebra/vc_vc.hpp"

#define __plugin algebra::vc
#define ALGEBRA_PLUGIN vc_array

namespace detray {

// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

// Define namespace(s)
namespace getter = algebra::getter;
namespace vector = algebra::vector;
namespace matrix = algebra::matrix;

// Define matrix operator
template <typename T>
using matrix_operator =
    matrix::actor<std::size_t, T, matrix::determinant::preset0<std::size_t, T>,
                  matrix::inverse::preset0<std::size_t, T>>;

}  // namespace detray
