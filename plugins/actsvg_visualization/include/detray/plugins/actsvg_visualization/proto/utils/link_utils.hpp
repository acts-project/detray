#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/surface_functors.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/definitions/math.hpp"


// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/core.hpp"

// System include(s)
#include <assert.h>
#include <cmath>

namespace detray::actsvg_visualization::proto::utils {

/// @brief Checks if the detray surface has a volume link.
template <typename detector_t>
inline auto has_link(const detray::surface<detector_t>& d_portal){
    const auto d_link_idx = d_portal.template visit_mask<get_link_functor>();
    return !is_invalid_value(d_link_idx);
}
/// @note expects that the detray surface has a volume link.
/// @returns the volume link of the detray surface.
template <typename detector_t>
inline auto get_link_volume(const detector_t& detector, const detray::surface<detector_t>& d_portal){
    assert(has_link(d_portal));
    const auto d_link_idx = d_portal.template visit_mask<get_link_functor>();
    return detector.volume_by_index(d_link_idx);
}

/// @brief Calculates the start and end point of the link.
/// @note expects that the detray surface has a volume link.
/// @returns (start, end).
template <typename detector_t>
inline auto link_points(const typename detector_t::geometry_context& context, const detector_t& detector, const detray::surface<detector_t>& d_portal, typename detector_t::point3 dir, const double link_length)
{
    assert(has_link(d_portal));

    // Calculating the start position:
    const auto start = d_portal.template visit_mask<utils::link_start_functor>(d_portal.transform(context));
    
    // Calculating the end position:
    const auto n = d_portal.normal(context, d_portal.global_to_local(context, start, dir));
    const auto volume_link = get_link_volume(detector, d_portal);
    const auto end = d_portal.template visit_mask<utils::link_end_functor>(detector, volume_link, start, n, link_length);

    return std::make_tuple(start, end);
}

}