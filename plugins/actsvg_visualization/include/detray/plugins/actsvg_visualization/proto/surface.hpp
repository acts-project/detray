#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/surface_functors.hpp"

namespace detray::actsvg_visualization::proto {

/// @brief Calculates the proto surface of a surface.
///
/// @param d_surface The detray surface.
/// @param context The context.
///
/// @note The transform is not taken into account for objects such as rings,
/// cylinders etc (not implemented yet).
///
/// @returns An actsvg proto surface representing the surface.
template <typename detector_t>
auto surface(const typename detector_t::geometry_context& context,
             const detray::surface<detector_t>& d_surface) {
    return d_surface.template visit_mask<utils::to_proto_surface_functor>(
        context, d_surface);
}

}  // namespace detray::actsvg_visualization::proto