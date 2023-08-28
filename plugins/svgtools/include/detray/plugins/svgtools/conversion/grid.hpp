/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"
#include "detray/grids/axis.hpp"
#include "detray/definitions/units.hpp"


// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <vector>
#include <algorithm>

namespace detray::svgtools::conversion {

template <typename point3_t>
inline std::string point_to_string(point3_t point) {

    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << "(" << point[0] << ", "
           << point[1] << ", " << point[2] << ")";
    return stream.str();
}

template <typename point3_t, typename detector_t>
inline auto information_section(
    const typename detector_t::geometry_context& context,
    const detray::surface<detector_t>& d_surface) {
    svgtools::meta::proto::information_section<point3_t> is;
    if (d_surface.is_portal()) {
        is._title = "Portal";
    } else {
        is._title = "Surface";
    }
    const auto position = d_surface.transform(context).translation();
    is._info = {"Idx: " + std::to_string(d_surface.index()),
                point_to_string(position)};
    is._position = svgtools::conversion::point<point3_t>(position);
    return is;
}

}  // namespace detray::svgtools::conversion