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

using point3 = std::array<actsvg::scalar, 3>;
using point3_container = std::vector<point3>;
using proto_surface = actsvg::proto::surface<point3_container>;
using proto_portal = actsvg::proto::portal<point3_container>;
using proto_link = proto_portal::link;

struct get_link {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index) const {
        const auto& m = mask_group[index];
        return m.volume_link();
    }
};

/// @returns returns the actsvg proto link from detray portal to volume.
template <typename detector_t>
proto_link create_p_link(const detray::surface<detector_t>& d_portal, const detray::detector_volume<detector_t>& d_volume, const typename detector_t::geometry_context& context){
    //actsvg::views::x_y view;
    proto_link p_link;
    const auto start = d_portal.center(context); //d_volume.center()

    

    const auto end = d_portal.center(context);//algebra::cmath::operator+(d_portal.center(context), d_portal.normal(context)); //+ d_portal.normal(context)
    //std::cout << typeid(d_volume.center()).name();
    //std::cout << typeid(d_portal.normal(context)).name();
    //std::cout << typeid(d_volume.center() + d_portal.normal(context)).name();
    p_link._start = convert_point<3>(start);
    p_link._end = convert_point<3>(end);//convert_point<3>(volume_position);
    auto p = end;
    std::cout<<"(" + std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " + std::to_string(p[3]) + ")";
    return p_link;
}

/// @returns An actsvg proto portal representing the portal.
/// @note detray portal is_portal() must be true.
template <typename detector_t>
proto_portal convert_portal(const detector_t& detector, const detray::surface<detector_t>& d_portal, const typename detector_t::geometry_context& context)
{
    assert(d_portal.is_portal());
    proto_portal p_portal;

    proto_surface p_surface = convert_surface(d_portal, context);
    p_portal._surface = p_surface;

    const auto d_link_idx = d_portal.template visit_mask<get_link>();

    // Check if there is a link
    if (d_link_idx != std::numeric_limits<decltype(d_link_idx)>::max())
    {
        const auto d_volume = detector.volume_by_index(d_link_idx);
        const auto p_link = create_p_link(d_portal, d_volume, context);
        p_portal._volume_links = std::vector{p_link};
    }

    return p_portal;
}
}