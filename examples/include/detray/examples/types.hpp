/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// This file contains the definitions of the linear algebra types to be used in
/// the examples, which in this case will be the std::array based plugin.
/// The includes for the different algebra plugins can be found in
/// @c detray/plugins/algebra/

// Detray core include(s).
#include "detray/plugins/algebra/array_definitions.hpp"

namespace detray::example {

using transform3 = algebra::array::transform3<detray::scalar>;
using point2 = algebra::array::point2<detray::scalar>;
using point3 = algebra::array::point3<detray::scalar>;
using vector2 = algebra::array::vector2<detray::scalar>;
using vector3 = algebra::array::vector3<detray::scalar>;

}  // namespace detray::example
