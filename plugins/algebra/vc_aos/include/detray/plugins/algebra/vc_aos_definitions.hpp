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

/// The plugin definition
template <algebra::concepts::scalar scalar_t>
using vc_aos = algebra::plugin::vc_aos<scalar_t>;

namespace getter {

using algebra::vc_aos::storage::block;
using algebra::vc_aos::storage::element;
using algebra::vc_aos::storage::set_block;
using algebra::vc_aos::storage::vector;

}  // namespace getter

namespace vector {

// Vc array specific
using algebra::vc_aos::math::cross;
using algebra::vc_aos::math::dot;
using algebra::vc_aos::math::eta;
using algebra::vc_aos::math::norm;
using algebra::vc_aos::math::normalize;

// No specific vectorized implementation needed
using algebra::vc_aos::math::perp;
using algebra::vc_aos::math::phi;
using algebra::vc_aos::math::theta;

}  // namespace vector

namespace matrix {

using algebra::vc_aos::math::identity;
using algebra::vc_aos::math::set_identity;
using algebra::vc_aos::math::set_zero;
using algebra::vc_aos::math::transpose;
using algebra::vc_aos::math::zero;

// Placeholder, until vectorization-friendly version is available
using algebra::vc_aos::math::determinant;
using algebra::vc_aos::math::inverse;

}  // namespace matrix

}  // namespace detray
