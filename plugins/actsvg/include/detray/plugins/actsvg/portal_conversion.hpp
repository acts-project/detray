#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/surface_conversion.hpp"

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
proto_link convert_link(const detray::surface<detector_t>& d_portal, const detray::detector_volume<detector_t>& d_volume, const typename detector_t::geometry_context& context){
    //actsvg::views::x_y view;
    //proto_link p_link;
    //const auto start = d_portal.center(context); //d_volume.center()

    
    //std::array point{0., 0., 0.};
    //const auto end = d_portal.normal(context, point);//algebra::cmath::operator+(d_portal.center(context), d_portal.normal(context)); //+ d_portal.normal(context)
    //std::cout << typeid(d_volume.center()).name();
    //std::cout << typeid(d_portal.normal(context)).name();
    //std::cout << typeid(d_volume.center() + d_portal.normal(context)).name();
    //p_link._start = convert_point<3>(std::array{0.,0.,0.});
    //p_link._end = convert_point<3>(std::array{0.,0.,100.});//convert_point<3>(volume_position);
    //auto p = start;
    //std::cout<<"(" + std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " + std::to_string(p[2]) + ") \n";

    //auto p2 = end;
    //std::cout<<"(" + std::to_string(p2[0]) + ", " + std::to_string(p2[1]) + ", " + std::to_string(p2[2]) + ")\n\n";

    const auto dir = typename detector_t::point3{};

    std::array point{0., 100., 0.};
    point = d_portal.global_to_local(context, point, dir);
    std::cout << "HERE:\n";
    const auto rim = d_portal.local_to_global(context, d_portal.closest_surface_point(point), dir);
    const auto n = d_portal.normal(context, point) * 100.;

    const auto start = rim;//std::array{10., 0., 0.};
    const auto end = n + start;
    std::cout << "start: " + point_to_string(start) + "\n";
    std::cout << "end: " + point_to_string(end) + "\n";
    std::cout << "rim: " + point_to_string(rim) + "\n";
    std::cout << "normal: " + point_to_string(n) + "\n";
    std::cout << "\n";
    return convert_link(start, end);
}

/// @returns An actsvg proto portal representing the portal.
/// @note detray portal is_portal() must be true.
template <typename detector_t>
proto_portal convert_portal(const detector_t& detector, const detray::surface<detector_t>& d_portal, const typename detector_t::geometry_context& context)
{
    assert(d_portal.is_portal());
    proto_portal p_portal;
    if (has_link(d_portal))
    {
        const auto d_volume = get_link_volume(detector, d_portal);
        const auto p_link = convert_link(d_portal, d_volume, context);
        p_portal._volume_links = {p_link};
    }
    p_portal._surface = convert_surface(d_portal, context);
    return p_portal;
}
}