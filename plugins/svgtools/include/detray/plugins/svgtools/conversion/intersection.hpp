/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/intersection.hpp"
#include "detray/plugins/svgtools/conversion/landmark.hpp"
#include "detray/plugins/svgtools/meta/proto/intersection_record.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools::conversion {

/// @returns The proto intersection of a detray intersection.
template <typename point3_t, typename detector_t, typename intersection_t>
inline auto intersection(
    const detector_t& detector, const dvector<intersection_t>& intersections,
    const typename detector_t::vector3& dir = {},
    const typename detector_t::geometry_context& gctx = {}) {
    using p_intersection_t = svgtools::meta::proto::intersection<point3_t>;
    p_intersection_t p_ir;

    std::vector<typename p_intersection_t::landmark_type> landmarks;
    for (const auto& intr : intersections) {
        const detray::surface surface{detector, intr.sf_desc};
        const auto position = surface.local_to_global(gctx, intr.local, dir);
        const auto p_lm = svgtools::conversion::landmark<point3_t>(position);
        landmarks.push_back(p_lm);
    }

    p_ir._landmarks = std::move(landmarks);
    return p_ir;
}

}  // namespace detray::svgtools::conversion
