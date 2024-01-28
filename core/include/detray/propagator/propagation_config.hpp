/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/navigation/navigation_config.hpp"
#include "detray/propagator/stepping_config.hpp"

namespace detray::propagation {

/// Configuration of the propagation
template <typename scalar_t>
struct config {
    navigation::config<scalar_t> navigation{};
    stepping::config<scalar_t> stepping{};
};

}  // namespace detray::propagation
