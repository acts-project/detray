/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/masks/annulus2D.hpp"
#include "detray/tools/generators.hpp"

// System include(s)
#include <cmath>
#include <iostream>

namespace detray::svgtools::utils {

template <typename mask_t>
inline auto local_vertices(const mask_t& mask) {
    return vertices<typename mask_t::point2_t, typename mask_t::point3_t>(mask,
                                                                          10);
}

inline auto local_vertices(const mask<line<>>& m) {
    using scalar_t = typename mask<line<>>::scalar_type;
    using point3_t = typename mask<line<>>::point3_t;
    const auto hz = m[line<>::e_half_z];
    scalar_t zero{0};
    point3_t p1{zero, zero, -hz};
    point3_t p2{zero, zero, hz};
    return dvector<point3_t>{p1, p2};
}

template <typename transform_t, typename mask_t>
inline auto global_vertices(const transform_t& transform, const mask_t& mask) {
    auto ret = local_vertices(mask);
    for (size_t i = 0; i < ret.size(); i++) {
        ret[i] = transform.point_to_global(ret[i]);
    }
    return ret;
}

}  // namespace detray::svgtools::utils