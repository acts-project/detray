#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/transform_conversion.hpp"
#include "detray/plugins/actsvg/proto_types.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

namespace {

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

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape An annulus2D.
template <typename detector_t, typename bounds_t>
auto convert_surface(const detray::surface<detector_t>& d_surface, const typename detector_t::geometry_context&, const detray::annulus2D<>& shape, const bounds_t& bounds)
{
    //Rotation for circular objects is currently not supported.
    assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});

    proto_surface p_surface;
    auto ri = static_cast<actsvg::scalar>(bounds[shape.e_min_r]);
    auto ro = static_cast<actsvg::scalar>(bounds[shape.e_max_r]);
    //p_surface._type = proto_surface::type::e_annulus;
    p_surface._type = proto_surface::type::e_disc;
    p_surface._radii = {ri, ro};
    //p_surface._zparameters = {nhz + hz, hz};
    return p_surface;
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape A cylinder2D.
template <typename detector_t, typename bounds_t, bool kRadialCheck, template <typename> class intersector_t>
auto convert_surface(const detray::surface<detector_t>& d_surface, const typename detector_t::geometry_context&, const detray::cylinder2D<kRadialCheck, intersector_t>& shape, const bounds_t& bounds)
{
    //Rotation for circular objects is currently not supported.
    assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});

    proto_surface p_surface;
    auto r = static_cast<actsvg::scalar>(bounds[shape.e_r]);
    auto nhz = static_cast<actsvg::scalar>(bounds[shape.e_n_half_z]);
    auto phz = static_cast<actsvg::scalar>(bounds[shape.e_p_half_z]);
    p_surface._type = proto_surface::type::e_cylinder;
    p_surface._radii = {0., r};
    auto hz = (phz-nhz)/2;
    p_surface._zparameters = {nhz + hz, hz};
    return p_surface;
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape A ring2D.
template <typename detector_t, typename bounds_t>
auto convert_surface(const detray::surface<detector_t>& d_surface, const typename detector_t::geometry_context& context, const detray::ring2D<>& shape, const bounds_t& bounds)
{
    //Rotation for circular objects is currently not supported.
    assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});

    proto_surface p_surface;
    auto ri = static_cast<actsvg::scalar>(bounds[shape.e_inner_r]);
    auto ro = static_cast<actsvg::scalar>(bounds[shape.e_outer_r]);
    auto center = convert_point<3>(d_surface.center(context));
    p_surface._type = proto_surface::type::e_disc;
    p_surface._radii = {ri, ro};
    p_surface._zparameters = {center[2], 0};
    return p_surface;
}

/// @brief Sets the proto surfaces type and bounds to be equivalent to the detray shape.
///
/// @param shape A polygon.
template <typename detector_t, typename bounds_t, typename polygon_t>
auto convert_surface(const detray::surface<detector_t>& d_surface, const typename detector_t::geometry_context& context, const polygon_t, const bounds_t&)
{
    proto_surface p_surface;
    p_surface._type = proto_surface::type::e_polygon;
    std::array dir{0.,0.,1.};
    set_vertices(p_surface, d_surface.global_vertices(context, dir));
    return p_surface;
}

/// @brief A functor to set the proto surfaces type and bounds to be equivalent to the mask.
struct convert_surface_functor {

    template <typename mask_group_t, typename index_t, typename detector_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index, const detray::surface<detector_t>& d_surface, const typename detector_t::geometry_context& context) const {
        const auto& m = mask_group[index];
        return convert_surface(d_surface, context, m.get_shape(), m.values());
    }
};

}
/// @brief Calculates the proto surface of a surface.
///
/// @param d_surface The detray surface.
/// @param context The context.
///
/// @note The transform is not taken into account for objects such as rings, cylinders etc (not implemented yet).
///
/// @returns An actsvg proto surface representing the surface.
template <typename detector_t>
auto convert_surface(
    const detray::surface<detector_t>& d_surface,
    const typename detector_t::geometry_context& context) {
    return d_surface.template visit_mask<convert_surface_functor>(d_surface, context);
}
}  // namespace detray::actsvg_visualization