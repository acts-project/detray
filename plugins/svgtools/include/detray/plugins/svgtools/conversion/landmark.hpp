/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/plugins/svgtools/conversion/point.hpp"
#include "detray/plugins/svgtools/meta/proto/landmark.hpp"
#include "detray/plugins/svgtools/utils/intersection_utils.hpp"

namespace detray::svgtools::conversion {

/// @returns The proto landmark of a detray intersection.
template <typename point3_t, typename detector_t>
inline auto landmark(
    const typename detector_t::geometry_context& context,
    const detector_t& detector,
    const detray::intersection2D<typename detector_t::surface_type,
                                 typename detector_t::transform3>&
        d_intersection) {
    const auto position =
        svgtools::utils::intersection_point(context, detector, d_intersection);
    return svgtools::conversion::landmark<point3_t>(position);
}

/// @returns The proto landmark of a detray point.
template <typename point3_t, typename d_point3_t>
inline auto landmark(d_point3_t& position) {
    using p_landmark_t = svgtools::meta::proto::landmark<point3_t>;
    p_landmark_t p_lm;
    p_lm._position = svgtools::conversion::point<point3_t>(position);
    return p_lm;
}

}  // namespace detray::svgtools::conversion
