/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/array_cmath.hpp"

namespace detray {

/// The plugin definition
template <typename scalar_t>
using array = algebra::plugin::array<scalar_t>;

using algebra::cmath::operator*;
using algebra::cmath::operator-;
using algebra::cmath::operator+;

namespace getter {

using algebra::cmath::storage::block;
using algebra::cmath::storage::element;
using algebra::cmath::storage::set_block;
using algebra::cmath::storage::vector;

}  // namespace getter

namespace vector {

// array specific implementations
using algebra::cmath::dot;
using algebra::cmath::normalize;

// generic implementations
using algebra::cmath::cross;
using algebra::cmath::eta;
using algebra::cmath::norm;
using algebra::cmath::perp;
using algebra::cmath::phi;
using algebra::cmath::theta;

}  // namespace vector

namespace matrix {

using algebra::cmath::identity;
using algebra::cmath::set_identity;
using algebra::cmath::set_zero;
using algebra::cmath::zero;

using algebra::cmath::determinant;
using algebra::cmath::inverse;
using algebra::cmath::transpose;

}  // namespace matrix

}  // namespace detray
