#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/surface.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/core.hpp"

namespace detray::actsvg_visualization::link {

namespace {
struct get_link {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index) const {
        const auto& m = mask_group[index];
        return m.volume_link();
    }
};

/// @brief Checks if the detray surface has a volume link.
template <typename detector_t>
auto has_link(const detray::surface<detector_t>& d_portal){
    const auto d_link_idx = d_portal.template visit_mask<get_link>();
    return d_link_idx != std::numeric_limits<decltype(d_link_idx)>::max();
}
/// @note expects that the detray surface has a volume link.
/// @returns the volume link of the detray surface.
template <typename detector_t>
auto get_link_volume(const detector_t& detector, const detray::surface<detector_t>& d_portal){
    const auto d_link_idx = d_portal.template visit_mask<get_link>();
    return detector.volume_by_index(d_link_idx);
}

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
auto closest_surface_point(const mask_t& mask, const transform_t& transform, const typename mask_t::point3_t& dir, const typename mask_t::point3_t& loc_calc_offset){
    const typename mask_t::local_frame_type frame{};
    const auto loc_center = bounding_box_center(mask);
    const auto point = frame.global_to_local(transform, loc_calc_offset + loc_center, dir);
    const auto loc_surface_point = mask.closest_surface_point(point);
    return frame.local_to_global(transform, loc_surface_point);
}

// Calculates the link starting location of the remaining shapes.
template <typename mask_t, typename transform_t>
auto link_start(const mask_t& mask, const transform_t& transform, const typename mask_t::point3_t& dir){
    return closest_surface_point(mask, transform, dir, {0.,0.,0.});
}

// Calculates the (optimal) link starting point for rings.
template <typename transform_t>
auto link_start(const detray::mask<detray::ring2D<>>& mask, const transform_t& transform, const typename detray::mask<detray::ring2D<>>::point3_t& dir){
    const auto shape = mask.get_shape();
    if (mask[shape.e_inner_r] == mask[shape.e_outer_r]){
        return closest_surface_point(mask, transform, dir, {0.,1.,0.});
    }
    const auto r = (mask[shape.e_inner_r] + mask[shape.e_outer_r])/2;
    const detray::mask<detray::ring2D<>> middle_circle{0u, r, r};
    return link_start(middle_circle, transform, dir);
}

// Calculates the (optimal) link starting point for cylinders.
template <typename transform_t, bool kRadialCheck, template <typename> class intersector_t>
auto link_start(const detray::mask<detray::cylinder2D<kRadialCheck, intersector_t>>& mask, const transform_t& transform, const typename detray::mask<detray::ring2D<>>::point3_t& dir){
return closest_surface_point(mask, transform, dir, {0.,1.,0.}); //TODO: Take into account z rotation.
}

struct link_start_functor {
    template <typename mask_group_t, typename index_t, typename transform_t, typename point_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index, const transform_t& transform, const point_t& dir) const {
        const auto& m = mask_group[index];
        return link_start(m, transform, dir);
    }
};

}

struct link_options{
    double length;
};

/// @returns The actsvg proto link from one point to another.
template <typename point3_t>
inline auto to_proto_link(const point3_t& start, const point3_t& end){
    conversion_types::link p_link;
    p_link._start = transform::convert_point<3>(start);
    p_link._end = transform::convert_point<3>(end);
    actsvg::style::color c{{255, 0, 0}, 0.75};
    auto marker = actsvg::style::marker({"<<"});
    marker._size = 1.2;
    p_link._stroke = actsvg::style::stroke(c);
    p_link._end_marker = marker;
    return p_link;
}

/// @returns a tuple with the link calculated using the surface normal. The first and second item is are the links in the negative and positive direction of the normal, respectively.
template <typename detector_t>
inline auto proto_links(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_portal, const link_options& l_options){
    typename detector_t::point3 dir{};
    const auto start = d_portal.template visit_mask<link_start_functor>(d_portal.transform(context), dir);
    const auto n = d_portal.normal(context, d_portal.global_to_local(context, start, dir));
    const auto pos_end = (n*l_options.length) + start;
    const auto neg_end = (n*-l_options.length) + start;
    return std::array{to_proto_link(start, neg_end), to_proto_link(start, pos_end)};
}

}