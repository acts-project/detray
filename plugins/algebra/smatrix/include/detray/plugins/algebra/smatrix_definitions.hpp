// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// TODO: Remove this when gcc fixes their false positives.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#endif

// Algebra-Plugins include
#include "algebra/smatrix_smatrix.hpp"

namespace detray {

/// The plugin definition
template <algebra::concepts::scalar scalar_t>
using smatrix = algebra::plugin::smatrix<scalar_t>;

namespace getter {

using algebra::smatrix::storage::block;
using algebra::smatrix::storage::element;
using algebra::smatrix::storage::set_block;
using algebra::smatrix::storage::vector;
}  // namespace getter

namespace vector {

using algebra::smatrix::math::cross;
using algebra::smatrix::math::dot;
using algebra::smatrix::math::eta;
using algebra::smatrix::math::norm;
using algebra::smatrix::math::normalize;
using algebra::smatrix::math::perp;
using algebra::smatrix::math::phi;
using algebra::smatrix::math::theta;

}  // namespace vector

namespace matrix {

using algebra::smatrix::math::determinant;
using algebra::smatrix::math::identity;
using algebra::smatrix::math::inverse;
using algebra::smatrix::math::set_identity;
using algebra::smatrix::math::set_zero;
using algebra::smatrix::math::transpose;
using algebra::smatrix::math::zero;

}  // namespace matrix

}  // namespace detray
