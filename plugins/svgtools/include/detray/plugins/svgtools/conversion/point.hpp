/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <algorithm>

namespace detray::svgtools::conversion {

template <typename ret_point_t, typename arg_point_t>
inline auto point(const arg_point_t& p) {
    ret_point_t ret{};
    const auto n = std::min(ret.size(), p.size());
    for (std::size_t i = 0; i < n; i++) {
        ret[i] = p[i];
    }
    return ret;
}

}  // namespace detray::svgtools::conversion