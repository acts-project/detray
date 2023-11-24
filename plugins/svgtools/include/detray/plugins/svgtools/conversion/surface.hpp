/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"

// Actsvg include(s)
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <algorithm>
#include <iterator>

namespace detray::svgtools::conversion {

/// @brief Sets the measures of the proto surface to be the same as the mask.
template <typename mask_t>
inline void set_measures(
    actsvg::proto::surface<std::vector<typename mask_t::point3_t>>& p_surface,
    const mask_t& m) {

    auto cast_scalar = [](const typename mask_t::scalar_type& v) {
        return static_cast<actsvg::scalar>(v);
    };
    std::transform(m.values().cbegin(), m.values().cend(),
                   std::back_inserter(p_surface._measures), cast_scalar);
}

/// @brief Sets the vertices of the proto surface to be the same as the mask.
template <typename transform_t, typename mask_t>
inline void set_vertices(
    actsvg::proto::surface<std::vector<typename transform_t::point3>>&
        p_surface,
    const transform_t& trf, const mask_t& m) {

    // Approximate any arcs in the mask shape with ten line segments
    auto vertices = m.vertices(10u);
    for (std::size_t i = 0; i < vertices.size(); i++) {
        p_surface._vertices.push_back(trf.point_to_global(vertices[i]));
    }
}

/// @brief Returns the proto surface for a shape.
/// @note For lines, the thickness is fixed and not determined by the cross
/// section.
template <typename transform_t, typename mask_t>
auto inline surface(const transform_t& transform, const mask_t& m) {

    using point3_t = typename transform_t::point3;
    using p_surface_t = actsvg::proto::surface<std::vector<point3_t>>;

    p_surface_t p_surface;
    p_surface._type = p_surface_t::type::e_polygon;
    set_measures(p_surface, m);
    set_vertices(p_surface, transform, m);

    return p_surface;
}

/// @brief Returns the proto surface for 2D cylinders.
template <typename transform_t, bool kRadialCheck,
          template <typename> class intersector_t>
auto inline surface(const transform_t& transform,
                    const mask<cylinder2D<kRadialCheck, intersector_t>>& m) {

    using shape_t = cylinder2D<kRadialCheck, intersector_t>;
    using point3_t = typename transform_t::point3;
    using p_surface_t = actsvg::proto::surface<std::vector<point3_t>>;

    p_surface_t p_surface;

    const auto r = static_cast<actsvg::scalar>(m[shape_t::e_r]);
    const auto nhz = static_cast<actsvg::scalar>(m[shape_t::e_n_half_z]);
    const auto phz = static_cast<actsvg::scalar>(m[shape_t::e_p_half_z]);
    const auto center = transform.translation();
    const auto hz = static_cast<actsvg::scalar>(0.5f * (phz - nhz) + center[2]);

    p_surface._type = p_surface_t::type::e_cylinder;
    p_surface._radii = {0.f, r};
    p_surface._zparameters = {nhz + hz, hz};
    p_surface._transform._tr = {static_cast<actsvg::scalar>(center[0]),
                                static_cast<actsvg::scalar>(center[1])};
    set_measures(p_surface, m);

    return p_surface;
}

/// @brief Returns the proto surface for 2D rings.
template <typename transform_t>
auto surface(const transform_t& transform, const mask<ring2D<>>& m) {

    using shape_t = ring2D<>;
    using point3_t = typename transform_t::point3;
    using p_surface_t = actsvg::proto::surface<std::vector<point3_t>>;

    p_surface_t p_surface;

    const auto ri = static_cast<actsvg::scalar>(m[shape_t::e_inner_r]);
    const auto ro = static_cast<actsvg::scalar>(m[shape_t::e_outer_r]);
    const auto center = transform.translation();
    const auto z = static_cast<actsvg::scalar>(center[2]);

    p_surface._type = p_surface_t::type::e_disc;
    p_surface._radii = {ri, ro};
    p_surface._zparameters = {z, z};
    p_surface._transform._tr = {static_cast<actsvg::scalar>(center[0]),
                                static_cast<actsvg::scalar>(center[1])};
    set_measures(p_surface, m);

    return p_surface;
}

/// @brief Returns the proto surface for 2D annuli.
template <typename transform_t>
auto inline surface(const transform_t& transform, const mask<annulus2D<>>& m) {

    using point3_t = typename transform_t::point3;
    using p_surface_t = actsvg::proto::surface<std::vector<point3_t>>;

    p_surface_t p_surface;

    p_surface._type = p_surface_t::type::e_annulus;
    set_measures(p_surface, m);
    set_vertices(p_surface, transform, m);

    return p_surface;
}

/// @brief Returns the proto surface for 2D rings.
template <typename transform_t, bool kSquareCrossSect,
          template <typename> class intersector_t>
auto surface(const transform_t& transform,
             const mask<line<kSquareCrossSect, intersector_t>>& m) {

    using shape_t = line<kSquareCrossSect, intersector_t>;
    using point3_t = typename transform_t::point3;
    using p_surface_t = actsvg::proto::surface<std::vector<point3_t>>;

    p_surface_t p_surface;

    // All line surfaces are drawn as a circles(straws) in xy-view
    const auto r{static_cast<actsvg::scalar>(m[shape_t::e_cross_section])};
    const auto hz{static_cast<actsvg::scalar>(m[shape_t::e_half_z])};
    const auto center = transform.translation();

    p_surface._type = p_surface_t::type::e_straw;
    p_surface._radii = {1.f, r};
    p_surface._zparameters = {-hz, hz};
    p_surface._transform._tr = {static_cast<actsvg::scalar>(center[0]),
                                static_cast<actsvg::scalar>(center[1])};
    set_measures(p_surface, m);

    return p_surface;
}

/// @brief Returns the proto surface for a shape.
/// @note For lines, the thickness is fixed and not determined by the cross
/// section.
template <typename point3_container_t, typename mask_t>
auto inline surface(const mask_t& m) {
    using transform_t = typename mask_t::local_frame_type::transform3_type;
    return detray::svgtools::conversion::surface(transform_t{}, m);
}

namespace {
/// @brief A functor to set the proto surfaces type and bounds to be equivalent
/// to the mask.
struct surface_converter {
    template <typename mask_group_t, typename index_t, typename transform_t>
    DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                       const index_t& index,
                                       const transform_t& transform) const {

        return svgtools::conversion::surface(transform, mask_group[index]);
    }
};
}  // namespace

/// @brief Calculates the proto surface of a surface.
///
/// @param d_surface The detray surface.
/// @param context The geometry context.
///
/// @note The transform is not taken into account for objects such as rings,
/// cylinders etc (not implemented yet).
///
/// @returns An actsvg proto surface representing the surface.
template <typename detector_t>
auto surface(const typename detector_t::geometry_context& context,
             const detray::surface<detector_t>& d_surface,
             const styling::surface_style& style =
                 styling::tableau_colorblind::surface_style_sensitive) {

    auto p_surface = d_surface.template visit_mask<surface_converter>(
        d_surface.transform(context));

    p_surface._name = "surface_" + std::to_string(d_surface.index());
    svgtools::styling::apply_style(p_surface, style);

    return p_surface;
}

}  // namespace detray::svgtools::conversion
