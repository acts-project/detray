#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/surface_functors.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/transform_utils.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/core.hpp"

namespace detray::actsvg_visualization::proto {

/// @returns The link calculated using the surface normal vector.
template <typename detector_t>
inline auto link(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_portal, const double link_length){
    typename detector_t::point3 dir{};
    const auto start = d_portal.template visit_mask<utils::link_start_functor>(d_portal.transform(context), dir);
    const auto n = d_portal.normal(context, d_portal.global_to_local(context, start, dir));
    const auto end = (n*link_length) + start;
    
    proto_link p_link;
    p_link._start = utils::convert_point<3>(start);
    p_link._end = utils::convert_point<3>(end);

    return p_link;
}
}