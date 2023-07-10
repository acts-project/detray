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

    template <class mask_t>
    inline proto_surface convert(const mask_t& m){
        proto_surface p_surface;

        const auto& shape = m.get_shape();

            constexpr auto is_annulus2D = std::is_same_v<decltype(m), const detray::mask<detray::annulus2D<>>&>;
            constexpr auto is_cuboid3D = std::is_same_v<decltype(m), const detray::mask<detray::cuboid3D<>>&>;
            constexpr auto is_cylinder2D = std::is_same_v<decltype(m), const detray::mask<detray::cylinder2D<>>&>;
            constexpr auto is_cylinder3D = std::is_same_v<decltype(m), const detray::mask<detray::cylinder3D>&>;
            constexpr auto is_line = std::is_same_v<decltype(m), const detray::mask<detray::line<>>&>;
            //constexpr auto is_masks = std::is_same_v<decltype(m), const detray::mask<detray::masks<>>&>;
            constexpr auto is_rectangle2D = std::is_same_v<decltype(m), const detray::mask<detray::rectangle2D<>>&>;
            constexpr auto is_ring2D = std::is_same_v<decltype(m), const detray::mask<detray::ring2D<>>&>;
            constexpr auto is_single3D = std::is_same_v<decltype(m), const detray::mask<detray::single3D<>>&>;
            constexpr auto is_trapezoid2D = std::is_same_v<decltype(m), const detray::mask<detray::trapezoid2D<>>&>;
            //constexpr auto is_unbounded = std::is_same_v<decltype(m), const detray::mask<detray::unbounded<>>&>;
            //constexpr auto is_unmasked = std::is_same_v<decltype(m), const detray::mask<detray::unmasked<>>>;

            if constexpr (is_annulus2D)
            {
                p_surface._type = proto_surface::type::e_annulus;
                auto ri = static_cast<scalar>(m[shape.e_min_r]);
                auto ro = static_cast<scalar>(m[shape.e_max_r]);

                p_surface._radii = {ri, ro};
            }

            else if constexpr (is_cylinder2D)
            {
                p_surface._type = proto_surface::type::e_cylinder;
                auto r = static_cast<scalar>(m[m.get_shape().e_r]);
                auto nhz = static_cast<scalar>(m[shape.e_n_half_z]);
                auto phz = static_cast<scalar>(m[shape.e_p_half_z]);
                
                p_surface._radii = {0., r};
                p_surface._zparameters = {-nhz, -phz};
            }

            else if constexpr (is_cylinder3D)
            {
                p_surface._type = proto_surface::type::e_cylinder;
                auto ri = static_cast<scalar>(m[shape.e_min_r]);
                auto ro = static_cast<scalar>(m[shape.e_max_r]);
                auto min_z = static_cast<scalar>(m[m.get_shape().e_min_z]);
                auto max_z = static_cast<scalar>(m[m.get_shape().e_max_z]);
                auto hz = static_cast<scalar>((max_z - min_z)/2.0);

                p_surface._radii = {ri, ro};
                p_surface._zparameters = {-hz, hz};
            }

            else if constexpr (is_line)
            {
                
            }

            else if constexpr (is_ring2D)
            {
                p_surface._type = proto_surface::type::e_disc;
                auto ri = static_cast<scalar>(m[shape.e_inner_r]);
                auto ro = static_cast<scalar>(m[shape.e_outer_r]);
                p_surface._radii = {ri, ro};
            }

            else if constexpr (is_rectangle2D || is_trapezoid2D || is_cuboid3D) {
                p_surface._type = proto_surface::type::e_polygon;

                auto detray_vertices = shape.template local_vertices<>(m.values());            
                point3_container actsvg_vertices;
                for (auto dv : detray_vertices)
                {
                    point3 av = {static_cast<scalar>(dv[0]), static_cast<scalar>(dv[1]), static_cast<scalar>(dv[2])};
                    actsvg_vertices.push_back(av);
                }

                p_surface._vertices = actsvg_vertices;
            }
            

            return p_surface;
    }

    template <typename algebra_t = __plugin::transform3<detray::scalar>>
    struct to_proto_surface {

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST inline proto::surface<point3_container> operator()(
            const mask_group_t& mask_group, 
            const index_t& index) const {

            const auto& m = mask_group[index];

            return convert(m);
        }
    };
    
}