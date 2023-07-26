#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/mask_conversion.hpp"
#include "detray/plugins/actsvg/transform_conversion.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

using point3 = std::array<actsvg::scalar, 3>;
using point3_container = std::vector<point3>;
using proto_surface = actsvg::proto::surface<point3_container>;


/// @brief A functor to set the proto surfaces type and bounds to be equivalent to the mask.
struct set_type_and_bounds_functor {

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline void operator()(
        const mask_group_t& mask_group, const index_t& index, proto_surface& p_surface) const {
        const auto& m = mask_group[index];
        set_type_and_bounds(p_surface, m.get_shape(), m.values());
    }
};

/// @brief Calculates the proto surface of a surface.
///
/// @param d_surface The detray surface.
/// @param context The context.
///
/// @note The transform is not taken into account for objects such as rings, cylinders etc (not implemented yet).
///
/// @returns An actsvg proto surface representing the surface.
template <typename detector_t>
proto_surface convert_surface(
    const detray::surface<detector_t>& d_surface,
    const typename detector_t::geometry_context& context) {
        proto_surface p_surface;
        d_surface.template visit_mask<set_type_and_bounds_functor>(p_surface);
        set_vertices(p_surface, d_surface.global_vertices(context));
        return p_surface;
}
}  // namespace detray::actsvg_visualization