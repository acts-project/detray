#include <vector>
#include <type_traits>

#include "detray/geometry/surface.hpp"
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

using namespace actsvg;

namespace surface_converter{

    using point3 = std::array<scalar, 3>;
    using point3_container = std::vector<point3>;
    using proto_surface = proto::surface<point3_container>;

    template <typename transform_t>
    inline auto convert_transform(transform_t d_transform){
        auto translation = d_transform.translation();
        auto euler_angles = rotation_matrix_to_euler_angles(d_transform.rotation);
        //TODO: skew and scale

        auto ret = actsvg::style::transform();
        ret._tr = {translation[0], translation[1]};
        ret._rot = {euler_angles[0], euler_angles[2], euler_angles[3]};

        return ret;
    }

    template <typename point3_t, typename matrix_t>
    inline point3_t rotation_matrix_to_euler_angles(matrix_t mat){
        float a = sqrt(mat[0][0] * mat[0][0] + mat[1][0] * mat[1][0]);
        // Checking if is singular.
        if (a < 1e-6)
            return {atan2(-mat[1][2], mat[1][1]), atan2(-mat[2][0], a), 0};

        return {atan2(mat[2][1], mat[2][2]), atan2(-mat[2][0], a), atan2(mat[1][0], mat[0][0])};
    }

    template <class mask_t>
    inline proto_surface convert_mask(const mask_t& m){
        proto_surface p_surface;

        const auto& shape = m.get_shape();

        constexpr auto is_annulus2D = std::is_same_v<decltype(m), const detray::mask<detray::annulus2D<>>&>;
        constexpr auto is_cylinder2D = std::is_same_v<decltype(m), const detray::mask<detray::cylinder2D<>>&>;
        constexpr auto is_rectangle2D = std::is_same_v<decltype(m), const detray::mask<detray::rectangle2D<>>&>;
        constexpr auto is_ring2D = std::is_same_v<decltype(m), const detray::mask<detray::ring2D<>>&>;
        constexpr auto is_trapezoid2D = std::is_same_v<decltype(m), const detray::mask<detray::trapezoid2D<>>&>;

        // Set bounds.
        if constexpr (is_annulus2D)
        {
            auto ri = static_cast<scalar>(m[shape.e_min_r]);
            auto ro = static_cast<scalar>(m[shape.e_max_r]);

            p_surface._type = proto_surface::type::e_annulus;
            p_surface._radii = {ri, ro};
        }

        else if constexpr (is_cylinder2D)
        {
            auto r = static_cast<scalar>(m[m.get_shape().e_r]);
            auto nhz = static_cast<scalar>(m[shape.e_n_half_z]);
            auto phz = static_cast<scalar>(m[shape.e_p_half_z]);
            
            p_surface._type = proto_surface::type::e_cylinder;
            p_surface._radii = {0., r};
            p_surface._zparameters = {-nhz, -phz};
        }

        else if constexpr (is_ring2D)
        {
            auto ri = static_cast<scalar>(m[shape.e_inner_r]);
            auto ro = static_cast<scalar>(m[shape.e_outer_r]);

            p_surface._type = proto_surface::type::e_disc;
            p_surface._radii = {ri, ro};
        }

        else if constexpr (is_rectangle2D)
        {
            p_surface._type = proto_surface::type::e_polygon;
        }

        else if constexpr (is_trapezoid2D)
        {
            p_surface._type = proto_surface::type::e_polygon;
        }

        // Set vertices.
        auto detray_vertices = shape.template local_vertices<>(m.values());            
        point3_container actsvg_vertices;
        for (auto dv : detray_vertices)
        {
            point3 av = {static_cast<scalar>(dv[0]), static_cast<scalar>(dv[1]), static_cast<scalar>(dv[2])};
            actsvg_vertices.push_back(av);
        }
        p_surface._vertices = actsvg_vertices;

        return p_surface;
    }

    template <typename algebra_t = __plugin::transform3<detray::scalar>>
    struct to_proto_surface {

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST inline proto::surface<point3_container> operator()(
            const mask_group_t& mask_group, 
            const index_t& index) const {
            const auto& m = mask_group[index];
            return convert_mask(m);
        }
    };
    
    template <typename surface_t, typename context_t>
    proto_surface convert_surface(surface_t d_surface, context_t context)
    {
        auto p_surface = d_surface.template visit_mask<to_proto_surface>();
        p_surface._transform = convert_transform(p_surface.transform(context));
        return p_surface;
    }
}