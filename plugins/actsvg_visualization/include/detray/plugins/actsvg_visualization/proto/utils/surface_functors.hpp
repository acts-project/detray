#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/surface_utils.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/transform_utils.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/core.hpp"

namespace detray::actsvg_visualization::proto::utils {

/// @returns functor to obtain the volume link.
struct get_link_functor {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index) const {
        const auto& m = mask_group[index];
        return m.volume_link();
    }
};

/// @returns functor to an optimal starting point for displaying the link.
struct link_start_functor {

    public:
    
    template <typename mask_group_t, typename index_t, typename transform_t, typename point_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index, const transform_t& transform, const point_t& dir) const {
        const auto& m = mask_group[index];
        return link_start(m, transform, dir);
    }

    private:

    // Calculates the link starting location of the remaining shapes.
    template <typename mask_t, typename transform_t>
    auto link_start(const mask_t& mask, const transform_t& transform, const typename mask_t::point3_t& dir) const {
        return nearest_point_from_center(mask, transform, dir, {0.,0.,0.});
    }

    // Calculates the (optimal) link starting point for rings.
    template <typename transform_t>
    auto link_start(const detray::mask<detray::ring2D<>>& mask, const transform_t& transform, const typename detray::mask<detray::ring2D<>>::point3_t& dir) const {
        const auto shape = mask.get_shape();
        if (mask[shape.e_inner_r] == mask[shape.e_outer_r]){
            return nearest_point_from_center(mask, transform, dir, {0.,1.,0.});
        }
        const auto r = (mask[shape.e_inner_r] + mask[shape.e_outer_r])/2;
        const detray::mask<detray::ring2D<>> middle_circle{0u, r, r};
        return link_start(middle_circle, transform, dir);
    }

    // Calculates the (optimal) link starting point for cylinders (2D).
    template <typename transform_t, bool kRadialCheck, template <typename> class intersector_t>
    auto link_start(const detray::mask<detray::cylinder2D<kRadialCheck, intersector_t>>& mask, const transform_t& transform, const typename detray::mask<detray::ring2D<>>::point3_t& dir) const {
    return nearest_point_from_center(mask, transform, dir, {0.,1.,0.}); //TODO: Take into account z rotation.
    }

};

/// @brief A functor to set the proto surfaces type and bounds to be equivalent to the mask.
struct to_proto_surface_functor {

    public:

    template <typename mask_group_t, typename index_t, typename detector_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index, const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface) const {
        const auto& m = mask_group[index];
        return to_proto_surface(context, d_surface, m.get_shape(), m.values());
    }

    private:

    // Returns the proto surface for remaining shapes.
    template <typename detector_t, typename bounds_t, typename polygon_t>
    auto to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const polygon_t, const bounds_t&) const 
    {
        proto::proto_surface p_surface;
        p_surface._type = proto::proto_surface::type::e_polygon;
        std::array dir{0.,0.,1.};
        set_vertices(p_surface, d_surface.global_vertices(context, dir));
        return p_surface;
    }

    // Returns the proto surface for 2D cylinders.
    template <typename detector_t, typename bounds_t, bool kRadialCheck, template <typename> class intersector_t>
    auto to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const detray::cylinder2D<kRadialCheck, intersector_t>& shape, const bounds_t& bounds) const 
    {
        //Rotation for circular objects is currently not supported.
        assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});
        //Only translation one z axis is supported.
        assert(d_surface.transform.translate(context).x() == 0 && d_surface.transform.translate(context).y() == 0);

        proto::proto_surface p_surface;
        auto r = static_cast<actsvg::scalar>(bounds[shape.e_r]);
        auto nhz = static_cast<actsvg::scalar>(bounds[shape.e_n_half_z]);
        auto phz = static_cast<actsvg::scalar>(bounds[shape.e_p_half_z]);
        p_surface._type = proto::proto_surface::type::e_cylinder;
        p_surface._radii = {static_cast<actsvg::scalar>(0), r};
        auto hz = (phz-nhz)/2 + static_cast<actsvg::scalar>(d_surface.center(context)[2]);
        p_surface._zparameters = {nhz + hz, hz};
        return p_surface;
    }

    // Returns the proto surface for 2D rings.
    template <typename detector_t, typename bounds_t>
    auto to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const detray::ring2D<>& shape, const bounds_t& bounds) const 
    {
        //Rotation for circular objects is currently not supported.
        assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});
        //Only translation one z axis is supported.
        assert(d_surface.transform.translate(context).x() == 0 && d_surface.transform.translate(context).y() == 0);

        proto::proto_surface p_surface;
        auto ri = static_cast<actsvg::scalar>(bounds[shape.e_inner_r]);
        auto ro = static_cast<actsvg::scalar>(bounds[shape.e_outer_r]);
        auto center = utils::convert_point<3>(d_surface.center(context));
        p_surface._type = proto::proto_surface::type::e_disc;
        p_surface._radii = {ri, ro};
        p_surface._zparameters = {center[2], static_cast<actsvg::scalar>(0)};
        return p_surface;
    }

    // Returns the proto surface for 2D annuluses.
    template <typename detector_t, typename bounds_t>
    auto to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const detray::annulus2D<>& shape, const bounds_t& bounds) const 
    {
        //Rotation for circular objects is currently not supported.
        assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});
        //Only translation one z axis is supported.
        assert(d_surface.transform.translate(context).x() == 0 && d_surface.transform.translate(context).y() == 0);

        proto::proto_surface p_surface;
        auto ri = static_cast<actsvg::scalar>(bounds[shape.e_min_r]);
        auto ro = static_cast<actsvg::scalar>(bounds[shape.e_max_r]);
        auto center = utils::convert_point<3>(d_surface.center(context));

        p_surface._type = proto_surface::type::e_annulus;
        p_surface._radii = {ri, ro};
        p_surface._zparameters = {center[2], static_cast<actsvg::scalar>(0)};

        std::array dir{0.,0.,1.};
        set_vertices(p_surface, d_surface.global_vertices(context, dir));

        return p_surface;
    }

    /// @brief Sets the vertices of the proto surfaces type to be equivalent to the detray shape.
    template <class container_t>
    void set_vertices(proto::proto_surface& p_surface, const container_t& vertices) const 
    {
        proto::point3_container actsvg_vertices;
        for (auto v : vertices) {
            actsvg_vertices.push_back(utils::convert_point<3>(v));
        }
        p_surface._vertices = actsvg_vertices;
    }
};

}