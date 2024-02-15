/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray core include(s).
#include "detray/definitions/detail/algebra.hpp"

namespace detray::test {

using transform3 = __plugin::transform3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;
using vector2 = __plugin::vector2<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;

#if DETRAY_ALGEBRA_ARRAY
static constexpr char filenames[] = "array-";
#elif DETRAY_ALGEBRA_EIGEN
static constexpr char filenames[] = "eigen-";
#elif DETRAY_ALGEBRA_SMATRIX
static constexpr char filenames[] = "smatrix-";
#elif DETRAY_ALGEBRA_VC
static constexpr char filenames[] = "vc-";
#endif

}  // namespace detray::test
