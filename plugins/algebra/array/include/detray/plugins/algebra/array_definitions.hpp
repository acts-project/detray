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
using point3D = transform3D::point3;
using vector3D = transform3D::vector3;
/// @}

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct array {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using transform3D = algebra::array::transform3<value_type>;
    using point2D = algebra::array::point2<value_type>;
    using point3D = algebra::array::point3<value_type>;
    using vector3D = algebra::array::vector3<value_type>;
};
/// @}

// Define namespace(s)
namespace matrix = algebra::matrix;

namespace vector {

using algebra::cmath::cross;
using algebra::cmath::dot;
using algebra::cmath::normalize;

}  // namespace vector

namespace getter {

using algebra::cmath::eta;
using algebra::cmath::norm;
using algebra::cmath::perp;
using algebra::cmath::phi;
using algebra::cmath::theta;

using algebra::cmath::element;

/// Function extracting a slice from a matrix
template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS,
          typename scalar_t>
ALGEBRA_HOST_DEVICE inline auto vector(
    const algebra::array::matrix_type<scalar_t, ROWS, COLS>& m, std::size_t row,
    std::size_t col) {

    return algebra::cmath::vector_getter<
        std::size_t, algebra::array::storage_type, scalar_t, SIZE>()(m, row,
                                                                     col);
}

}  // namespace getter

// Define matrix/vector operator
template <typename scalar_t>
using standard_matrix_operator =
    matrix::actor<scalar_t, matrix::determinant::preset0<scalar_t>,
                  matrix::inverse::preset0<scalar_t>>;

}  // namespace detray
