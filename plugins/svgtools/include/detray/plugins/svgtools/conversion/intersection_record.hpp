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
#include "detray/plugins/svgtools/utils/intersection_utils.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools::conversion {

/// @returns The proto intersection record from a list of points
template <typename point3_t, typename container_t>
inline auto intersection_record(const container_t& points) {
    using p_intersection_record_t =
        svgtools::meta::proto::intersection_record<point3_t>;
    p_intersection_record_t p_ir;
    for (const auto& point : points) {
        const auto p_lm = svgtools::conversion::landmark<point3_t>(point);
        p_ir._landmarks.push_back(p_lm);
    }
    return p_ir;
}

/// @returns The proto intersection record of a detray intersection record.
template <typename point3_t, typename detector_t>
inline auto intersection_record(
    const typename detector_t::geometry_context& context,
    const detector_t& detector,
    const std::vector<
        std::pair<detray::dindex,
                  detray::intersection2D<typename detector_t::surface_type,
                                         typename detector_t::transform3>>>&
        intersection_record) {
    std::vector<typename detector_t::point3> points;
    for (const auto& record : intersection_record) {
        const auto point = svgtools::utils::intersection_point(
            context, detector, record.second);
        points.push_back(point);
    }
    return svgtools::conversion::intersection_record<point3_t>(points);
}

}  // namespace detray::svgtools::conversion