#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/surface_conversion.hpp"
#include "detray/definitions/units.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/core.hpp"

// Algebra include(s)
#include "algebra/math/cmath.hpp"

// System include(s)
#include <assert.h>
#include <vector>
#include <type_traits>
#include <iostream>

namespace detray::actsvg_visualization {

namespace {
struct get_link {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index) const {
        const auto& m = mask_group[index];
        return m.volume_link();
    }
};

template <typename detector_t>
auto has_link(const detray::surface<detector_t>& d_portal){
    const auto d_link_idx = d_portal.template visit_mask<get_link>();
    return d_link_idx != std::numeric_limits<decltype(d_link_idx)>::max();
}

template <typename detector_t>
auto get_link_volume(const detector_t& detector, const detray::surface<detector_t>& d_portal){
    const auto d_link_idx = d_portal.template visit_mask<get_link>();
    return detector.volume_by_index(d_link_idx);
}

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

template <typename mask_t, typename transform_t>
auto closest_surface_point(const mask_t& mask, const transform_t& transform, const typename mask_t::point3_t& dir, const typename mask_t::point3_t& loc_calc_offset){
    const typename mask_t::local_frame_type frame{};
    const auto loc_center = bounding_box_center(mask);
    const auto point = frame.global_to_local(transform, loc_calc_offset + loc_center, dir);
    const auto loc_surface_point = mask.closest_surface_point(point);
    return frame.local_to_global(transform, loc_surface_point);
}

template <typename mask_t, typename transform_t, typename view_t>
auto link_start(const mask_t& mask, const transform_t& transform, const typename mask_t::point3_t& dir, const view_t& view){
    if (std::is_same_v<decltype(view), const actsvg::views::x_y&>){
        return closest_surface_point(mask, transform, dir, {0.,1.,0.});
    }
    if (std::is_same_v<decltype(view), const actsvg::views::z_r&>){
        return closest_surface_point(mask, transform, dir, {0.,1.,0.});
    }
    return closest_surface_point(mask, transform, dir, {0.,0.,0.});
}

template <typename transform_t, typename view_t>
auto link_start(const detray::mask<detray::ring2D<>>& mask, const transform_t& transform, const typename detray::mask<detray::ring2D<>>::point3_t& dir, const view_t& view){
    const auto shape = mask.get_shape();
    if (mask[shape.e_inner_r] == mask[shape.e_outer_r]){
        return closest_surface_point(mask, transform, dir, {0.,1.,0.});
    }
    const auto r = (mask[shape.e_inner_r] + mask[shape.e_outer_r])/2;
    const detray::mask<detray::ring2D<>> middle_circle{0u, r, r};
    return link_start(middle_circle, transform, dir, view);
}

struct link_start_functor {
    template <typename mask_group_t, typename index_t, typename transform_t, typename point_t, typename view_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index, const transform_t& transform, const point_t& dir, const view_t& view) const {
        const auto& m = mask_group[index];
        return link_start(m, transform, dir, view);
    }
};

template <typename detector_t, typename view_t>
auto link_start(const detray::surface<detector_t>& d_portal, const typename detector_t::geometry_context& context, const typename detector_t::point3& dir, const view_t& view){
    return d_portal.template visit_mask<link_start_functor>(d_portal.transform(context), dir, view);
}

}

/// @returns returns the actsvg proto link from detray portal to volume.
template <typename point3_t>
proto_link convert_link(const point3_t& start, const point3_t& end){
    proto_link p_link;
    p_link._start = convert_point<3>(start);
    p_link._end = convert_point<3>(end);
    actsvg::style::color c{{255, 0, 0}, 0.9};
    p_link._stroke = actsvg::style::stroke(c);
    return p_link;
}
template <typename point3_t>
std::string point_to_string(point3_t point){
    return "(" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ", " + std::to_string(point[2]) + ")";
}

/// @returns returns the actsvg proto link from detray portal to volume.
template <typename detector_t>
proto_link convert_link(const detray::surface<detector_t>& d_portal, const typename detector_t::geometry_context& context){
    actsvg::views::x_y view;
    typename detector_t::point3 dir{};
    const auto start = link_start(d_portal, context, dir, view);
    const auto n = d_portal.normal(context, d_portal.global_to_local(context, start, dir));
    const auto end = (n*10.) + start;
    std::cout << "HERE: \n";
    std::cout << "start: " + point_to_string(start) + "\n";
    std::cout << "end: " + point_to_string(end) + "\n";
    std::cout << "normal: " + point_to_string(n) + "\n";
    std::cout << "\n";
    return convert_link(start, end);
}

}