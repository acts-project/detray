#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/conversion_types.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/transform_utils.hpp"
#include "detray/definitions/units.hpp"

// System include(s)
#include <optional>

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
    
    template <typename mask_group_t, typename index_t, typename transform_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index, const transform_t& transform) const {
        const auto& m = mask_group[index];
        return link_start(m, transform);
    }

    private:

    // Calculates the link starting location of the remaining shapes.
    template <typename mask_t, typename transform_t>
    auto inline link_start(const mask_t& mask, const transform_t& transform) const {
        return mask.local_min_bounds_center(transform);
    }

    // Calculates the (optimal) link starting point for rings.
    template <typename transform_t>
    auto inline link_start(const detray::mask<detray::ring2D<>>& mask, const transform_t& transform) const {
        using mask_t = typename detray::mask<detray::ring2D<>>;
        using shape_t = typename mask_t::shape;
        using point3_t = typename mask_t::point3_t;
        using scalar_t = typename mask_t::scalar_type;

        const scalar_t r{(mask[shape_t::e_inner_r] + mask[shape_t::e_outer_r])/2};
        const scalar_t phi{detray::constant<scalar_t>::pi_2};
        const scalar_t z{0};

        // Polar coordinate system.
        const typename mask_t::local_frame_type frame{};
        
        const auto true_center = mask.local_min_bounds_center(transform);
        const auto rel_point = frame.local_to_global(transform, point3_t{r, phi, z}) - transform.translation();
        return rel_point + true_center;
    }

    // Calculates the (optimal) link starting point for annuluses.
    template <typename transform_t>
    auto inline link_start(const detray::mask<detray::annulus2D<>>& mask, const transform_t& transform) const {
        using mask_t = typename detray::mask<detray::annulus2D<>>;
        using shape_t = typename mask_t::shape;
        using point3_t = typename mask_t::point3_t;
        using scalar_t = typename mask_t::scalar_type;

        const scalar_t r{(mask[shape_t::e_min_r] + mask[shape_t::e_max_r])/2};
        const scalar_t phi{mask[shape_t::e_average_phi]};
        const scalar_t z{0};
        
        // Polar coordinate system.
        const typename mask_t::local_frame_type frame{};

        const auto true_center = mask.local_min_bounds_center(transform);
        const auto rel_point = frame.local_to_global(transform, point3_t{r, phi, z}) - transform.translation();
        return rel_point + true_center;
    }

    // Calculates the (optimal) link starting point for cylinders (2D).
    template <typename transform_t, bool kRadialCheck, template <typename> class intersector_t>
    auto inline link_start(const detray::mask<detray::cylinder2D<kRadialCheck, intersector_t>>& mask, const transform_t& transform) const {
        using mask_t = typename detray::mask<detray::cylinder2D<kRadialCheck, intersector_t>>;
        using shape_t = typename mask_t::shape;
        using point3_t = typename mask_t::point3_t;
        using scalar_t = typename mask_t::scalar_type;

        const scalar_t r{mask[shape_t::e_r]};
        const scalar_t phi{detray::constant<scalar_t>::pi_2};
        const scalar_t z{0};

        // Cylindrical coordinate system.
        const typename mask_t::local_frame_type frame{};

        const auto true_center = mask.local_min_bounds_center(transform);        
        const auto rel_point = frame.local_to_global(transform, point3_t{r*phi, z, r}) - transform.translation();
        return rel_point + true_center;
    }

    // Calculates the (optimal) link starting point for cylinders (3D).
    template <typename transform_t>
    auto inline link_start(const detray::mask<detray::cylinder3D>& mask, const transform_t& transform) const {
        using mask_t = typename detray::mask<detray::cylinder3D>;
        using shape_t = typename mask_t::shape;
        using point3_t = typename mask_t::point3_t;
        using scalar_t = typename mask_t::scalar_type;

        const scalar_t r{(mask[shape_t::e_min_r] + mask[shape_t::e_max_r])/2};
        const scalar_t phi{(mask[shape_t::e_max_phi] + mask[shape_t::e_max_phi])/2};
        const scalar_t z{0};

        // Cylindrical coordinate system.
        const typename mask_t::local_frame_type frame{};

        const auto true_center = mask.local_min_bounds_center(transform);        
        const auto rel_point = frame.local_to_global(transform, point3_t{r*phi, z, r}) - transform.translation();
        return rel_point + true_center;
    }

};

/// @brief Functor to calculate the outermost radius a shape.
/// If the shape is not defined by a radius, then null option is returned.
struct outer_radius_functor
{

    public:

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline std::optional<detray::scalar> operator()(
    const mask_group_t& mask_group, const index_t& index) const {
        const auto& m = mask_group[index];
        return outer_radius(m);
    }

    private:

    // The remaining shapes do not have an outer radius.
    template <typename mask_t>
    std::optional<detray::scalar> inline outer_radius(const mask_t& /*mask*/) const {
        return std::nullopt;
    }

    // Calculates the outer radius for rings.
    auto inline outer_radius(const detray::mask<detray::ring2D<>>& mask) const {
        return std::optional<detray::scalar>(mask[mask.get_shape().e_outer_r]);
    }

    // Calculates the outer radius for annuluses.
    auto inline outer_radius(const detray::mask<detray::annulus2D<>>& mask) const {
        return std::optional<detray::scalar>(mask[mask.get_shape().e_max_r]);
    }

    // Calculates the outer radius for cylinders (2D).
    template <bool kRadialCheck, template <typename> class intersector_t>
    auto inline outer_radius(const detray::mask<detray::cylinder2D<kRadialCheck, intersector_t>>& mask) const {
        return std::optional<detray::scalar>(mask[mask.get_shape().e_r]);
    }

    // Calculates the outer radius for cylinders (3D).
    auto inline outer_radius(const detray::mask<detray::cylinder3D>& mask) const {
        return std::optional<detray::scalar>(mask[mask.get_shape().e_max_r]);
    }
    
};

struct link_end_functor {

    public:
    
    template <typename mask_group_t, typename index_t, typename detector_t, typename point3_t, typename scalar_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index, const detector_t& detector, const detray::detector_volume<detector_t>& volume_link, const point3_t& surface_point, const point3_t& surface_normal, const scalar_t& link_length) const {
        const auto& m = mask_group[index];
        return link_dir(m, detector, volume_link, surface_point, surface_normal) * link_length + surface_point;
    }

    private: 

    /// @brief Calculates the direction of the link for remaining shapes.
    template <typename detector_t, typename mask_t, typename point3_t>
    inline auto link_dir(const mask_t& /*mask*/, const detector_t& /*detector*/, const detray::detector_volume<detector_t>& volume_link, const point3_t& surface_point, const point3_t& surface_normal) const
    {
        const auto dir = volume_link.center() - surface_point;
        const auto dot_prod = algebra::cmath::dot(dir, surface_normal);
        typename detector_t::scalar_type sgn{0};
        if (dot_prod > 0){
            sgn = 1;
        }
        if (dot_prod < 0){
            sgn = -1;
        }
        return sgn * surface_normal;
    }

    /// @brief Calculates the direction of the link for cylinders (2D)
    template <typename detector_t, bool kRadialCheck, template <typename> class intersector_t, typename point3_t>
    inline auto link_dir(const detray::mask<detray::cylinder2D<kRadialCheck, intersector_t>>& mask, const detector_t& detector, const detray::detector_volume<detector_t>& volume_link, const point3_t& /*surface_point*/, const point3_t& surface_normal) const
    {
        for (const auto& desc : volume_link.surface_lookup())
        {
            const detray::surface surface{detector, desc};
            if (surface.is_portal()){
                if (auto radius = surface.template visit_mask<outer_radius_functor>()){
                    if (*radius > mask[mask.get_shape().e_r]){
                        return surface_normal;
                    }
                }
            } 
        }
        return typename detector_t::scalar_type{-1} * surface_normal;
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

    /// @brief Returns the proto surface for remaining shapes.
    template <typename detector_t, typename bounds_t, typename polygon_t>
    auto inline to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const polygon_t, const bounds_t&) const 
    {
        proto::proto_surface p_surface;
        p_surface._type = proto::proto_surface::type::e_polygon;
        std::array dir{0.,0.,1.};
        set_vertices(p_surface, d_surface.global_vertices(context, dir));
        return p_surface;
    }

    /// @brief Returns the proto surface for 2D cylinders.
    template <typename detector_t, typename bounds_t, bool kRadialCheck, template <typename> class intersector_t>
    auto inline to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const detray::cylinder2D<kRadialCheck, intersector_t>& shape, const bounds_t& bounds) const 
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

    /// @brief Returns the proto surface for 2D rings.
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

    /// @brief Returns the proto surface for 2D annuluses.
    template <typename detector_t, typename bounds_t>
    auto inline to_proto_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface, const detray::annulus2D<>& shape, const bounds_t& bounds) const 
    {
        //Rotation for circular objects is currently not supported.
        assert(d_surface.transform.rotation(context) == typename detector_t::transform3{});
        
        //Only translation on z axis is supported.
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
    void inline set_vertices(proto::proto_surface& p_surface, const container_t& vertices) const 
    {
        proto::point3_container actsvg_vertices;
        for (auto v : vertices) {
            actsvg_vertices.push_back(utils::convert_point<3>(v));
        }
        p_surface._vertices = actsvg_vertices;
    }
};

}