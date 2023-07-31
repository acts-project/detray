// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/actsvg/portal_conversion.hpp"
#include "detray/plugins/actsvg/svg_conversion.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/core/defs.hpp"
#include "actsvg/display/geometry.hpp"
#include "actsvg/meta.hpp"

// Algebra include(s)
#include "algebra/math/cmath.hpp"

// System include(s)
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

using namespace actsvg;

template <typename point3_t>
std::string point_to_string(point3_t point){
    return "(" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ", " + std::to_string(point[2]) + ")";
}

template <typename detector_t>
auto bounding_box_center(const detray::surface<detector_t>& d_portal, const typename detector_t::geometry_context& context){
    const auto local_min_bounds = d_portal.local_min_bounds();
    const auto center = d_portal.center(context);
    const typename detector_t::point3 bounding_box_center{(local_min_bounds[3] + local_min_bounds[0])/2, (local_min_bounds[4] + local_min_bounds[1])/2, (local_min_bounds[5] + local_min_bounds[2])/2};
    return bounding_box_center + center;
}

int main(int, char**) {

    using toy_detector_t = detray::detector<detray::toy_metadata<>>;

    using point3 = std::array<actsvg::scalar, 3>;
    using point3_container = std::vector<point3>;
    using proto_portal = actsvg::proto::portal<point3_container>;

    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    toy_detector_t::geometry_context context{};

    views::z_r view;

    actsvg::style::color c{{255, 0, 0}, 0.9};

    // Draw x-y-axis.
    style::stroke stroke_black = style::stroke();
    auto axis =
        draw::x_y_axes("xy", {-250, 250}, {-250, 250}, stroke_black, "x", "y");
    std::vector<actsvg::svg::object> objs;
    objs.push_back(axis);
    for (size_t i = 0; i < det.portals().size(); i++){ 

    const auto portal_desc = det.portals()[i];
    const auto d_portal = detray::surface{det, portal_desc};

    const auto dir = toy_detector_t::point3{};


    toy_detector_t::point3 point{0., 1., 0};
    const auto center = bounding_box_center(d_portal, context);
    point = d_portal.global_to_local(context, point + center, dir);
    const auto rim = d_portal.local_to_global(context, d_portal.closest_surface_point(point), dir);
    const auto n = d_portal.normal(context, point) * 10.;

    const auto start = rim;//std::array{10., 0., 0.};
    const auto end = n + start;
    auto p_link = detray::actsvg_visualization::convert_link(d_portal, context);
    p_link._stroke = style::stroke(c);
    proto_portal p_portal;
    const auto svg = display::portal_link("link", p_portal, p_link, view);
    const auto portal_svg = display::surface("portal", detray::actsvg_visualization::convert_surface(d_portal, context), view);
    objs.push_back(svg);
    objs.push_back(portal_svg);
    std::cout << "start: " + point_to_string(start) + "\n";
    std::cout << "end: " + point_to_string(end) + "\n";
    std::cout << "center: " + point_to_string(center) + "\n";
    std::cout << "rim: " + point_to_string(rim) + "\n";
    std::cout << "point: " + point_to_string(point) + "\n";
    std::cout << "normal: " + point_to_string(n) + "\n";


    }
    detray::actsvg_visualization::write_svg("link.svg", objs);
}