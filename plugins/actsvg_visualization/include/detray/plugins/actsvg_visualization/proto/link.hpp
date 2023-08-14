/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/link_utils.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/surface_functors.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/transform_utils.hpp"

namespace detray::actsvg_visualization::proto {

/// @returns The link calculated using the surface normal vector.
template <typename detector_t>
inline auto link(const typename detector_t::geometry_context& context,
                 const detector_t& detector,
                 const detray::surface<detector_t>& d_portal) {
    typename detector_t::point3 dir{};

    // Length of link arrow is hardcoded to 3.
    constexpr double link_length = 3.;

    const auto [start, end] =
        proto::utils::link_points(context, detector, d_portal, dir, link_length);

    proto_link p_link;
    p_link._start = utils::convert_point<3>(start);
    p_link._end = utils::convert_point<3>(end);

    return p_link;
}

}  // namespace detray::actsvg_visualization::proto