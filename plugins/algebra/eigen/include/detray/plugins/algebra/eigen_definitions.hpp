/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/eigen_eigen.hpp"

namespace detray {

/// The plugin definition
template <typename scalar_t>
using eigen = algebra::plugin::eigen<scalar_t>;

namespace getter {

using algebra::eigen::storage::block;
using algebra::eigen::storage::element;
using algebra::eigen::storage::set_block;
using algebra::eigen::storage::vector;

}  // namespace getter

namespace vector {

using algebra::eigen::math::cross;
using algebra::eigen::math::dot;
using algebra::eigen::math::eta;
using algebra::eigen::math::norm;
using algebra::eigen::math::normalize;

using algebra::eigen::math::perp;
using algebra::eigen::math::phi;
using algebra::eigen::math::theta;

}  // namespace vector

namespace matrix {

using algebra::eigen::math::determinant;
using algebra::eigen::math::identity;
using algebra::eigen::math::inverse;
using algebra::eigen::math::set_identity;
using algebra::eigen::math::set_zero;
using algebra::eigen::math::transpose;
using algebra::eigen::math::zero;

}  // namespace matrix

}  // namespace detray
