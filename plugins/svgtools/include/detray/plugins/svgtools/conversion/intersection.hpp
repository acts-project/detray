/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/plugins/svgtools/conversion/landmark.hpp"
#include "detray/plugins/svgtools/meta/proto/intersection.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools::conversion {

/// @returns The proto intersection of a detray intersection.
template <typename detector_t, typename intersection_t>
inline auto intersection(const detector_t& detector,
                         const dvector<intersection_t>& intersections,
                         const typename detector_t::vector3_type& dir = {},
                         const typename detector_t::geometry_context& gctx = {},
                         const styling::landmark_style& style =
                             styling::svg_default::intersection_style) {

    using point3_t = typename detector_t::point3_type;
    using p_intersection_t = svgtools::meta::proto::intersection<point3_t>;
    p_intersection_t p_ir;

    for (const auto& intr : intersections) {
        const detray::surface surface{detector, intr.sf_desc};
        const auto position = surface.local_to_global(gctx, intr.local, dir);
        const auto p_lm = svgtools::conversion::landmark(position, style);
        p_ir._landmarks.push_back(p_lm);
    }

    svgtools::styling::apply_style(p_ir, style);

    return p_ir;
}

}  // namespace detray::svgtools::conversion
