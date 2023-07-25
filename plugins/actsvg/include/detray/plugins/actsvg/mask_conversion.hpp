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

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape An annulus2D.
template <typename bounds_t>
void set_type_and_bounds(proto_surface& p_surface, const detray::annulus2D<>& shape, const bounds_t& bounds)
{
    auto ri = static_cast<actsvg::scalar>(bounds[shape.e_min_r]);
    auto ro = static_cast<actsvg::scalar>(bounds[shape.e_max_r]);
    //p_surface._type = proto_surface::type::e_annulus;
    p_surface._type = proto_surface::type::e_disc;
    p_surface._radii = {ri, ro};
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape A cylinder2D.
template <typename bounds_t, bool kRadialCheck, template <typename> class intersector_t>
void set_type_and_bounds(proto_surface& p_surface, const detray::cylinder2D<kRadialCheck, intersector_t>& shape, const bounds_t& bounds)
{
        auto r = static_cast<actsvg::scalar>(bounds[shape.e_r]);
        auto nhz = static_cast<actsvg::scalar>(bounds[shape.e_n_half_z]);
        auto phz = static_cast<actsvg::scalar>(bounds[shape.e_p_half_z]);
        p_surface._type = proto_surface::type::e_cylinder;
        p_surface._radii = {0., r};
        p_surface._zparameters = {-nhz, -phz};
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape A cylinder2D.
template <typename bounds_t>
void set_type_and_bounds(proto_surface& p_surface, const detray::ring2D<>& shape, const bounds_t& bounds)
{
        auto ri = static_cast<actsvg::scalar>(bounds[shape.e_inner_r]);
        auto ro = static_cast<actsvg::scalar>(bounds[shape.e_outer_r]);
        p_surface._type = proto_surface::type::e_disc;
        p_surface._radii = {ri, ro};
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape A polygon.
template <typename bounds_t, typename polygon_t>
void set_type_and_bounds(proto_surface& p_surface, const polygon_t, const bounds_t&)
{
    p_surface._type = proto_surface::type::e_polygon;
}

/// @brief Sets the vertices of the proto surfaces type to be equivalent to the detray shape.
template <class container_t>
void set_vertices(proto_surface& p_surface, const container_t& vertices)
{
    point3_container actsvg_vertices;
    for (auto v : vertices) {
        actsvg_vertices.push_back(convert_point<3>(v));
    }
    p_surface._vertices = actsvg_vertices;
}

/// @brief Converts a detray mask to an actsvg proto surface.
///
/// @param mask The detray mask.
///
/// @returns An actsvg proto surface with the parameters of the mask.
template <class mask_t>
proto_surface convert_mask(const mask_t& mask) {
    proto_surface p_surface;
    set_type_and_bounds(p_surface, mask.get_shape(), mask.values());
    set_vertices(p_surface, mask.local_vertices());
    return p_surface;
}
}  // namespace detray::actsvg_visualization