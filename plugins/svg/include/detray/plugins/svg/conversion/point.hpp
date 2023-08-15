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

namespace detray::svg::conversion {

template <typename ret_point_t, typename arg_point_t>
inline auto point(const arg_point_t& p)
{
    ret_point_t ret{};
    const auto n = std::min(ret.size(), p.size());
    for (size_t i = 0; i < n; i++){
        ret[i] = p[i];
    }
    return ret;
}

}  // namespace detray::actsvg_visualization::proto::utils