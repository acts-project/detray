#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/core.hpp"

namespace detray::actsvg_visualization::proto::utils {

/// @returns The center point in local cartesian coordinates of the local min bounds of the mask.
template <typename mask_t>
auto bounding_box_center(const mask_t& mask){
    const auto b = mask.local_min_bounds();
    const auto s = b.get_shape();
    return typename mask_t::point3_t{
        (b[s.e_max_x] + b[s.e_min_x])/2,
        (b[s.e_max_y] + b[s.e_min_y])/2,
        (b[s.e_max_z] + b[s.e_min_z])/2
    };
}

/// @brief Calculates the closest point of the surface.
/// @note The point used to find the closest point is the sum of loc_calc_offset and the center of the shape's min bounds bounding box.
template <typename mask_t, typename transform_t>
auto nearest_point_from_center(const mask_t& mask, const transform_t& transform, const typename mask_t::point3_t& dir, const typename mask_t::point3_t& loc_calc_offset){
    const typename mask_t::local_frame_type frame{};
    const auto loc_center = bounding_box_center(mask);
    const auto point = frame.global_to_local(transform, loc_calc_offset + loc_center, dir);
    const auto loc_surface_point = mask.nearest_point(point);
    return frame.local_to_global(transform, loc_surface_point);
}

}