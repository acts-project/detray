#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/transform.hpp"
#include "detray/plugins/actsvg_visualization/link.hpp"
#include "detray/plugins/actsvg_visualization/conversion_types.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization::surface {

namespace {

/// @brief Sets the vertices of the proto surfaces type to be equivalent to the detray shape.
template <class container_t>
void set_vertices(conversion_types::surface& p_surface, const container_t& vertices)
{
    conversion_types::point3_container actsvg_vertices;
    for (auto v : vertices) {
        actsvg_vertices.push_back(transform::convert_point<3>(v));
    }
    p_surface._vertices = actsvg_vertices;
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape An annulus2D.
template <typename detector_t, typename bounds_t>
auto convert_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const detray::annulus2D<>& shape, const bounds_t& bounds)
{
    //Rotation for circular objects is currently not supported.
    assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});
    //Only translation one z axis is supported.
    assert(d_surface.transform.translate(context).x() == 0 && d_surface.transform.translate(context).y() == 0);

    conversion_types::surface p_surface;
    auto ri = static_cast<actsvg::scalar>(bounds[shape.e_min_r]);
    auto ro = static_cast<actsvg::scalar>(bounds[shape.e_max_r]);
    //p_surface._type = proto_surface::type::e_annulus;
    p_surface._type = conversion_types::surface::type::e_disc;
    p_surface._radii = {ri, ro};
    //p_surface._zparameters = {nhz + hz, hz};
    return p_surface;
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape A cylinder2D.
template <typename detector_t, typename bounds_t, bool kRadialCheck, template <typename> class intersector_t>
auto to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const detray::cylinder2D<kRadialCheck, intersector_t>& shape, const bounds_t& bounds)
{
    //Rotation for circular objects is currently not supported.
    assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});
    //Only translation one z axis is supported.
    assert(d_surface.transform.translate(context).x() == 0 && d_surface.transform.translate(context).y() == 0);

    conversion_types::surface p_surface;
    auto r = static_cast<actsvg::scalar>(bounds[shape.e_r]);
    auto nhz = static_cast<actsvg::scalar>(bounds[shape.e_n_half_z]);
    auto phz = static_cast<actsvg::scalar>(bounds[shape.e_p_half_z]);
    p_surface._type = conversion_types::surface::type::e_cylinder;
    p_surface._radii = {0., r};
    auto hz = (phz-nhz)/2;
    p_surface._zparameters = {nhz + hz, hz};
    return p_surface;
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape A ring2D.
template <typename detector_t, typename bounds_t>
auto to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const detray::ring2D<>& shape, const bounds_t& bounds)
{
    //Rotation for circular objects is currently not supported.
    assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});
    //Only translation one z axis is supported.
    assert(d_surface.transform.translate(context).x() == 0 && d_surface.transform.translate(context).y() == 0);

    conversion_types::surface p_surface;
    auto ri = static_cast<actsvg::scalar>(bounds[shape.e_inner_r]);
    auto ro = static_cast<actsvg::scalar>(bounds[shape.e_outer_r]);
    auto center = transform::convert_point<3>(d_surface.center(context));
    p_surface._type = conversion_types::surface::type::e_disc;
    p_surface._radii = {ri, ro};
    p_surface._zparameters = {center[2], 0};
    return p_surface;
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape A polygon.
template <typename detector_t, typename bounds_t, typename polygon_t>
auto to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const polygon_t, const bounds_t&)
{
    conversion_types::surface p_surface;
    p_surface._type = conversion_types::surface::type::e_polygon;
    std::array dir{0.,0.,1.};
    set_vertices(p_surface, d_surface.global_vertices(context, dir));
    return p_surface;
}

/// @brief A functor to set the proto surfaces type and bounds to be equivalent to the mask.
struct to_proto_surface_functor {

    template <typename mask_group_t, typename index_t, typename detector_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index, const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface) const {
        const auto& m = mask_group[index];
        return to_proto_surface(context, d_surface, m.get_shape(), m.values());
    }
};

}

struct surface_options
{

};

struct portal_options
{
    surface_options s_options;
    link::link_options l_options;
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
auto to_proto_surface(
    const typename detector_t::geometry_context& context,
    const detray::surface<detector_t>& d_surface,
    const surface_options& s_options
    ) {
    return d_surface.template visit_mask<to_proto_surface_functor>(context, d_surface);
}

/// @returns An actsvg proto portal representing the portal.
/// @note detray portal is_portal() must be true.
template <typename detector_t>
auto to_proto_portal(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_portal, const portal_options& p_options)
{
    assert(d_portal.is_portal());
    conversion_types::portal p_portal;
    if (link::has_link(d_portal))
    {
        const auto [l1, l2] = link::proto_links(context, d_portal, p_options.l_options);
        p_portal._volume_links = {l1, l2};
    }
    p_portal._surface = to_proto_surface(context, d_portal, p_options.s_options);
    return p_portal;
}

template <typename detector_t, typename view_t>
auto to_svg(const typename detector_t::geometry_context& context, const view_t& view, const detray::surface<detector_t>& d_surface, const surface_options& s_options, const portal_options& p_options, const std::string& name)
{
    if (d_surface.is_portal())
    {
        auto p_portal = to_proto_portal(context, d_surface, p_options);
        return actsvg::display::portal(name, p_portal, view);
    }
    auto p_surface = to_proto_surface(context, d_surface, s_options);
    return actsvg::display::surface(name, p_surface, view);
}

}  // namespace detray::actsvg_visualization