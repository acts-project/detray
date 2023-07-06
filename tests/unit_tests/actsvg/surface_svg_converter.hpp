#include <vector>
#include <type_traits>

#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

using namespace actsvg;

template <typename algebra_t = __plugin::transform3<detray::scalar>>
struct to_proto_surface {

    using point3 = std::array<scalar, 3>;
    using point3_container = std::vector<point3>;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline proto::surface<point3_container> operator()(
        const mask_group_t& mask_group, 
        const index_t& index) const {

        const auto& m = mask_group[index];
        proto::surface<point3_container> proto_surface;

        constexpr auto is_ring2D = std::is_same_v<decltype(m), const detray::mask<detray::ring2D<>>&>;
        constexpr auto is_cylinder2D = std::is_same_v<decltype(m), const detray::mask<detray::cylinder2D<>>&>;
        
        if constexpr (is_ring2D)
        {
            proto_surface._type = proto::surface<point3_container>::type::e_disc;
            proto_surface._radii = {
                static_cast<scalar>(m[m.get_shape().e_inner_r]), 
                static_cast<scalar>(m[m.get_shape().e_outer_r])
                };
        }

        else if constexpr (is_cylinder2D)
        {
            proto_surface._type = proto::surface<point3_container>::type::e_cylinder;
            proto_surface._radii = {
                0, 
                static_cast<scalar>(m[m.get_shape().e_r])
                };
        }

        else if constexpr (std::is_same_v<decltype(m), const detray::mask<detray::rectangle2D<>>&>)
        {
            detray::mask<detray::rectangle2D<>> u;
            auto vertices2 = u.get_shape().local_vertices<>(m.values());
            point3_container vertices3;
            for (auto p2 : vertices2)
            {
                point3 p3 = {static_cast<scalar>(p2[0]), static_cast<scalar>(p2[1]), 0};
                vertices3.push_back(p3);
            }

            proto_surface._vertices = vertices3;
            //auto c = m.get_shape().template vertices<algebra_t>(m.values());
        }
        
        else {
            proto_surface._type = proto::surface<point3_container>::type::e_polygon;
            //point3_container vertices = m.vertices();
            //proto_surface._vertices = vertices;
        }

        return proto_surface;        
    }
};