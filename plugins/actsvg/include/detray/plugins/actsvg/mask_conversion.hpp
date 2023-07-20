#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
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

template <class mask_t>
void set_proto_surface_type_and_bounds(proto_surface& p_surface, const mask_t& mask)
{
    const auto& shape = mask.get_shape();

    using shape_t = typename mask_t::shape;
    constexpr auto is_annulus2D = std::is_same_v<shape_t, detray::annulus2D<>>;
    constexpr auto is_cylinder2D =
        std::is_same_v<shape_t, detray::cylinder2D<>> || std::is_same_v<shape_t, detray::cylinder2D<false, detray::cylinder_portal_intersector>>;
    constexpr auto is_rectangle2D =
        std::is_same_v<shape_t, detray::rectangle2D<>>;
    constexpr auto is_ring2D = std::is_same_v<shape_t, detray::ring2D<>>;
    constexpr auto is_trapezoid2D =
        std::is_same_v<shape_t, detray::trapezoid2D<>>;

    // Set bounds.
    if constexpr (is_annulus2D) {
        auto ri = static_cast<actsvg::scalar>(mask[shape.e_min_r]);
        auto ro = static_cast<actsvg::scalar>(mask[shape.e_max_r]);

        //p_surface._type = proto_surface::type::e_annulus;
        p_surface._type = proto_surface::type::e_disc;
        p_surface._radii = {ri, ro};
    }

    else if constexpr (is_cylinder2D) {
        auto r = static_cast<actsvg::scalar>(mask[shape.e_r]);
        auto nhz = static_cast<actsvg::scalar>(mask[shape.e_n_half_z]);
        auto phz = static_cast<actsvg::scalar>(mask[shape.e_p_half_z]);

        p_surface._type = proto_surface::type::e_cylinder;
        p_surface._radii = {0., r};
        p_surface._zparameters = {-nhz, -phz};
    }

    else if constexpr (is_ring2D) {
        auto ri = static_cast<actsvg::scalar>(mask[shape.e_inner_r]);
        auto ro = static_cast<actsvg::scalar>(mask[shape.e_outer_r]);

        p_surface._type = proto_surface::type::e_disc;
        p_surface._radii = {ri, ro};
    }

    else if constexpr (is_rectangle2D) {
        p_surface._type = proto_surface::type::e_polygon;
    }

    else if constexpr (is_trapezoid2D) {
        p_surface._type = proto_surface::type::e_polygon;
    }
}

template <class container_t>
void set_proto_surface_vertices(proto_surface& p_surface, const container_t& vertices)
{
    point3_container actsvg_vertices;
    for (auto dv : detray_vertices) {
        actsvg_vertices.push_back(convert_point<3>(dv));
    }
    p_surface._vertices = actsvg_vertices;

    return p_surface;
}

/// @brief Converts a detray mask to an actsvg proto surface.
///
/// @param mask The detray mask.
///
/// @returns An actsvg proto surface with the parameters of the mask.
template <class mask_t>
inline proto_surface convert_mask(const mask_t& mask) {
    proto_surface p_surface;
    set_proto_surface_type_and_bounds(p_surface, mask);
    set_proto_surface_vertices(p_surface, mask.local_vertices());
    return p_surface;
}
}  // namespace detray::actsvg_visualization