/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/array.hpp"

namespace detray {

/// The plugin definition
template <algebra::concepts::scalar scalar_t>
using array = algebra::plugin::array<scalar_t>;

using algebra::array::operator*;
using algebra::array::operator-;
using algebra::array::operator+;

namespace getter {

using algebra::array::storage::block;
using algebra::array::storage::element;
using algebra::array::storage::set_block;
using algebra::array::storage::vector;

}  // namespace getter

namespace vector {

// array specific implementations
using algebra::array::dot;
using algebra::array::normalize;

// generic implementations
using algebra::array::cross;
using algebra::array::eta;
using algebra::array::norm;
using algebra::array::perp;
using algebra::array::phi;
using algebra::array::theta;

}  // namespace vector

namespace matrix {

using algebra::array::identity;
using algebra::array::set_identity;
using algebra::array::set_zero;
using algebra::array::zero;

using algebra::array::determinant;
using algebra::array::inverse;
using algebra::array::transpose;

}  // namespace matrix

}  // namespace detray
