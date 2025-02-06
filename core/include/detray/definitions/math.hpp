// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

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
