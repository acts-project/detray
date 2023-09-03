/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray core include(s).
#include "detray/definitions/algebra.hpp"

namespace detray::test {

using scalar = detray::scalar;
using transform3 = __plugin::transform3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;
using vector2 = __plugin::vector2<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;

#if DETRAY_ALGEBRA_ARRAY

// The std::array based algebra plugin is always available in the tests
template <typename T = test::scalar>
using algebra_t = detray::cmath<T>;
static constexpr char filenames[] = "array-";

#elif DETRAY_ALGEBRA_EIGEN

template <typename T = test::scalar>
using algebra_t = detray::eigen<T>;
static constexpr char filenames[] = "eigen-";

#elif DETRAY_ALGEBRA_SMATRIX

template <typename T = test::scalar>
using algebra_t = detray::smatrix<T>;
static constexpr char filenames[] = "smatrix-";

#elif DETRAY_ALGEBRA_VC

template <typename T = test::scalar>
using algebra_t = detray::vc<T>;
static constexpr char filenames[] = "vc-";
#endif

}  // namespace detray::test
