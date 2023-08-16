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

/// @returns The proto intersection record of a detray intersection record.
template <typename point3_t, typename detector_t>
inline auto intersection_record(const typename detector_t::geometry_context& context, const detector_t& detector, const std::vector<std::pair<detray::dindex, detray::intersection2D<typename detector_t::surface_type, typename detector_t::transform3>>>& intersection_record){
    using p_intersection_record_t = svgtools::meta::proto::intersection_record<point3_t>;
    p_intersection_record_t p_ir;
    std::vector<typename p_intersection_record_t::landmark_type> landmarks;
    for (const auto& intersection : intersection_record){
        const auto p_lm = svgtools::conversion::landmark<point3_t>(context, detector, intersection.second);
        landmarks.push_back(p_lm);
    }
    p_ir._landmarks = landmarks;
    return p_ir;
}

}