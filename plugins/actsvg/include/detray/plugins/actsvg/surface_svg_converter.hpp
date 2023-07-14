#include <math.h>

#include <type_traits>
#include <vector>

// Project include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"

namespace detray::actsvg_visualization {

using point3 = std::array<actsvg::scalar, 3>;
using point3_container = std::vector<point3>;
using proto_surface = actsvg::proto::surface<point3_container>;

template <typename matrix_t>
inline std::array<detray::scalar, 3> rotation_matrix_to_euler_angles(
    const matrix_t& mat) {
    float a = std::sqrt(mat[0][0] * mat[0][0] + mat[1][0] * mat[1][0]);
    // Checking if it is singular.
    if (a < 1e-6)
        return {std::atan2(-mat[1][2], mat[1][1]), std::atan2(-mat[2][0], a),
                0};

    return {std::atan2(mat[2][1], mat[2][2]), std::atan2(-mat[2][0], a),
            std::atan2(mat[1][0], mat[0][0])};
}

template <typename transform_t>
inline auto convert_transform(const transform_t& d_transform) {
    auto translation = d_transform.translation();
    auto euler_angles =
        rotation_matrix_to_euler_angles<>(d_transform.rotation());

    // TODO: skew and scale

    auto ret = actsvg::style::transform();
    constexpr auto rad_to_deg = 180.0 / 3.14;
    ret._tr = {static_cast<actsvg::scalar>(translation[0]),
               static_cast<actsvg::scalar>(translation[1])};
    ret._rot = {static_cast<actsvg::scalar>(euler_angles[2] * rad_to_deg),
                static_cast<actsvg::scalar>(euler_angles[1] * rad_to_deg),
                static_cast<actsvg::scalar>(euler_angles[0] * rad_to_deg)};

    return ret;
}

template <class mask_t>
inline proto_surface convert_mask(const mask_t& m) {
    proto_surface p_surface;

    const auto& shape = m.get_shape();

    using shape_t = typename mask_t::shape;
    constexpr auto is_annulus2D = std::is_same_v<shape_t, detray::annulus2D<>>;
    constexpr auto is_cylinder2D =
        std::is_same_v<shape_t, detray::cylinder2D<>>;
    constexpr auto is_rectangle2D =
        std::is_same_v<shape_t, detray::rectangle2D<>>;
    constexpr auto is_ring2D = std::is_same_v<shape_t, detray::ring2D<>>;
    constexpr auto is_trapezoid2D =
        std::is_same_v<shape_t, detray::trapezoid2D<>>;

    // Set bounds.
    if constexpr (is_annulus2D) {
        auto ri = static_cast<actsvg::scalar>(m[shape.e_min_r]);
        auto ro = static_cast<actsvg::scalar>(m[shape.e_max_r]);

        p_surface._type = proto_surface::type::e_annulus;
        p_surface._radii = {ri, ro};
    }

    else if constexpr (is_cylinder2D) {
        auto r = static_cast<actsvg::scalar>(m[shape.e_r]);
        auto nhz = static_cast<actsvg::scalar>(m[shape.e_n_half_z]);
        auto phz = static_cast<actsvg::scalar>(m[shape.e_p_half_z]);

        p_surface._type = proto_surface::type::e_cylinder;
        p_surface._radii = {0., r};
        p_surface._zparameters = {-nhz, -phz};
    }

    else if constexpr (is_ring2D) {
        auto ri = static_cast<actsvg::scalar>(m[shape.e_inner_r]);
        auto ro = static_cast<actsvg::scalar>(m[shape.e_outer_r]);

        p_surface._type = proto_surface::type::e_disc;
        p_surface._radii = {ri, ro};
    }

    else if constexpr (is_rectangle2D) {
        p_surface._type = proto_surface::type::e_polygon;
    }

    else if constexpr (is_trapezoid2D) {
        p_surface._type = proto_surface::type::e_polygon;
    }

    // Set vertices.
    auto detray_vertices = m.local_vertices();
    point3_container actsvg_vertices;
    for (auto dv : detray_vertices) {
        point3 av = {static_cast<actsvg::scalar>(dv[0]),
                     static_cast<actsvg::scalar>(dv[1]),
                     static_cast<actsvg::scalar>(dv[2])};
        actsvg_vertices.push_back(av);
    }
    p_surface._vertices = actsvg_vertices;

    return p_surface;
}

struct to_proto_surface {

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline actsvg::proto::surface<point3_container> operator()(
        const mask_group_t& mask_group, const index_t& index) const {
        const auto& m = mask_group[index];
        return convert_mask(m);
    }
};

template <typename detector_t>
proto_surface convert_surface(
    const detray::surface<detector_t>& d_surface,
    const typename detector_t::geometry_context& context) {
    auto p_surface = d_surface.template visit_mask<to_proto_surface>();
    p_surface._transform = convert_transform(d_surface.transform(context));
    return p_surface;
}
}  // namespace detray::actsvg_visualization