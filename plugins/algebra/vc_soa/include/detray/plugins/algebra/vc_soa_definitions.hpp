/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/vc_soa.hpp"

namespace detray {

/// The plugin definition
template <typename scalar_t>
using vc_soa = algebra::plugin::vc_soa<scalar_t>;

// Pull in additional arithmetic operators for the algebra types
using algebra::storage::operator*;
using algebra::storage::operator/;
using algebra::storage::operator-;
using algebra::storage::operator+;

namespace detail {

// Pull in SoA overloads for boolean masks
using namespace ::algebra::boolean;

}  // namespace detail

namespace math {

// Import the overload between single value and simd math functions
using namespace ::algebra::math;

}  // namespace math

namespace getter {

using algebra::vc_soa::storage::block;
using algebra::vc_soa::storage::element;
using algebra::vc_soa::storage::set_block;
using algebra::vc_soa::storage::vector;

}  // namespace getter

namespace vector {

using algebra::vc_soa::math::cross;
using algebra::vc_soa::math::dot;
using algebra::vc_soa::math::eta;
using algebra::vc_soa::math::norm;
using algebra::vc_soa::math::normalize;
using algebra::vc_soa::math::perp;
using algebra::vc_soa::math::phi;
using algebra::vc_soa::math::theta;

}  // namespace vector

namespace matrix {

using algebra::vc_soa::math::determinant;
using algebra::vc_soa::math::identity;
using algebra::vc_soa::math::inverse;
using algebra::vc_soa::math::set_identity;
using algebra::vc_soa::math::set_zero;
using algebra::vc_soa::math::transpose;
using algebra::vc_soa::math::zero;

}  // namespace matrix

}  // namespace detray
