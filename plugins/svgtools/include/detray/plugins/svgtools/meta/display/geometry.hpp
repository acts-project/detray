/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project inlude(s)
#include "detray/plugins/svgtools/meta/proto/intersection_record.hpp"
#include "detray/plugins/svgtools/meta/proto/landmark.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <string>
#include <vector>

namespace detray::svgtools::meta::display {

/// @brief Converts a proto landmark to a SVG object.
template <typename point3_t, typename view_t>
inline auto landmark(const std::string& id,
                     const svgtools::meta::proto::landmark<point3_t>& lm,
                     const view_t& view) {
    const auto point_view = view(std::vector{lm._position})[0];
    return actsvg::draw::marker(id, point_view, lm._marker);
}

/// @brief Converts a proto intersection record to a SVG object.
template <typename point3_t, typename view_t>
inline auto intersection_record(
    const std::string& id,
    const svgtools::meta::proto::intersection_record<point3_t>& ir,
    const view_t& view) {
    actsvg::svg::object ret{._tag = "g", ._id = id};
    for (size_t index = 0; index < ir._landmarks.size(); index++) {
        const auto lm = ir._landmarks[index];
        const auto svg = svgtools::meta::display::landmark(
            id + "_intersection" + std::to_string(index), lm, view);
        ret.add_object(svg);
    }
    return ret;
}

}  // namespace detray::svgtools::meta::display