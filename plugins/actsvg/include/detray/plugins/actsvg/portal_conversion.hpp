#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/surface_conversion.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

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

template <typename link_t>
proto_link convert_link(link_t d_link)
{
    proto_link p_link;

    actsvg::style::stroke stroke_black = actsvg::style::stroke();

    //p_link._start = convert_point<3>(/* start point of d_link */);
    //p_link._end = convert_point<3>(/* end point of d_link */);
    p_link._stroke = stroke_black;
    return p_link;
}

template <typename detector_t>
proto_portal convert_portal(const detector_t& detector, const detray::surface<detector_t>& d_portal, const typename detector_t::geometry_context& context)
{
    assert(d_portal.is_portal());
    proto_portal p_portal;
    proto_surface p_surface = convert_surface(d_portal, context);
    p_portal._surface = p_surface;
    const auto d_link_idx = d_portal.template visit_mask<get_link>();
    const auto d_link_volume = detector.volume_by_index(d_link_idx);



    std::cout << typeid(d_link_idx).name();
    return p_portal;

    //p_portal._surface = convert_surface(d_portal, context);
    //p_portal._volume_links = std::transform(/* portal links begin*/, /* portal links end */, p_portal._volume_links.begin(), p_portal._volume_links.end(), convert_link)
    /* for each d_link in d_portal convert it and include it in p_portal volume links. */
}
}