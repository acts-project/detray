/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Not used currently in visualization implementation i.e. detray-actsvg
// conversion.

#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"

// System include(s)
#include <math.h>

#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization::proto::utils {

/// @brief Calculates the detray point3 as an actsvg point.
///
/// @param d_point The detray point3.
///
/// @returns An actsvg point3.
template <std::size_t dim, typename point_t>
inline std::array<actsvg::scalar, dim> convert_point(const point_t& d_point) {
    std::array<actsvg::scalar, dim> ret;
    for (std::size_t i = 0; i < ret.size(); i++) {
        ret[i] = static_cast<actsvg::scalar>(d_point[i]);
    }
    return ret;
}
}  // namespace detray::actsvg_visualization::proto::utils