/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-plugins include(s)
#include "detray/definitions/algebra.hpp"

// SYCL include(s).
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

namespace detray {

namespace math {

/// Namespace to pick up math functions from
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
using namespace ::sycl;
#else
using namespace ::algebra::math;
#endif  // SYCL

}  // namespace math

namespace detail {

using math::copysign;
using math::signbit;

}  // namespace detail

}  // namespace detray
