// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <array>
#include <string>


#include "detray/intersection/detail/trajectories.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/particle_gun.hpp"
#include "detray/plugins/svgtools/conversion/point.hpp"

using point2 = std::array<actsvg::scalar, 2>;

template <typename point3_t>
inline auto point_to_view(const point3_t& point, const actsvg::views::x_y& /*view*/)
{
    return point;
}

template <typename point3_t>
inline auto point_to_view(const point3_t& point, const actsvg::views::z_r& /*view*/)
{
    return point3_t{point[2], algebra::getter::perp(point)};
}

template <typename detector_t, typename intersection_t, typename view_t>
actsvg::svg::object svg_intersection(const std::string& identification, const typename detector_t::geometry_context& context, const detector_t& detector, const intersection_t& intersection, const view_t& view)
{
    const typename detector_t::point3 dir{};
    const detray::surface surface{detector, intersection.sf_desc};
    const auto point_view = point_to_view(surface.local_to_global(context, intersection.local, dir), view);
    const auto actsvg_point_view = detray::svgtools::conversion::point<point2>(point_view);
    return actsvg::draw::circle(identification, actsvg_point_view, 1., actsvg::style::fill({{0, 255, 0}, 5.}));
}

template <typename detector_t, typename intersection_t, typename view_t>
actsvg::svg::object svg_trace(const std::string& identification, const typename detector_t::geometry_context& context, const detector_t& detector, const std::vector<std::pair<detray::dindex, intersection_t>>& intersection_record, const view_t& view)
{
    actsvg::svg::object ret{._tag = "g", ._id = identification};
    for (size_t index = 0; index < intersection_record.size(); index++){

        const auto svg = svg_intersection(identification + "_intersection" + std::to_string(index), context, detector, intersection_record[index].second, view);
        ret.add_object(svg);
    }
    return ret;
}

int main(int, char**) {

    // Axes.
    const auto axes =
        actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250}, actsvg::style::stroke(), "axis1", "axis2");

    // Creating the view.
    const actsvg::views::z_r view{};
    // Creating the detector and geomentry context.
    using detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    detector_t::geometry_context context{};

    // Svg for the detector.
    const detray::svgtools::illustrator il{det, names};
    const auto det_svg = il.draw_detector("detector", context, view);
    
    // Creating the rays.
    using intersection_t = detray::intersection2D<typename detector_t::surface_type,
                                          typename detector_t::transform3>;
    using transform3_t = typename detector_t::transform3;

    unsigned int theta_steps{10u};
    unsigned int phi_steps{10u};
    const typename detector_t::point3 ori{0.f, 0.f, 100.f};


    size_t index = 0;
    // Iterate through uniformly distributed momentum directions with ray
    for (const auto test_ray :
         detray::uniform_track_generator<detray::detail::ray<transform3_t>>(
             theta_steps, phi_steps, ori)) {

        // Record all intersections and objects along the ray
        const auto intersection_record =
            detray::particle_gun::shoot_particle(det, test_ray);

        const std::string name = "trace" + std::to_string(index);
        const auto svg_tra = svg_trace(name, context, det, intersection_record, view);

        detray::svgtools::write_svg(name + ".svg", {det_svg, svg_tra});

        index++;
    }


}