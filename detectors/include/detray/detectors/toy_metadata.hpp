/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/detectors/odd_metadata.hpp"

namespace detray {

/// Defines the data types needed for the toy detector (same as ODD)
template <concepts::algebra algebra_t>
using toy_metadata = odd_metadata<algebra_t>;

}  // namespace detray
