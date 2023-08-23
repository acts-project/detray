/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"
#include "detray/grids/axis.hpp"
#include "detray/definitions/units.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <array>
#include <string>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <exception>


std::string label_to_string(const detray::n_axis::label& label) {
    switch (label) {
        case detray::n_axis::label::e_x:
            return "e_x";
        case detray::n_axis::label::e_y:
            return "e_y";
        case detray::n_axis::label::e_z:
            return "e_z";
        /*case detray::n_axis::label::e_r:
            return "e_r";
        case detray::n_axis::label::e_phi:
            return "e_phi";
        case detray::n_axis::label::e_rphi:
            return "e_rphi";
        case detray::n_axis::label::e_cyl_z:
            return "e_cyl_z";*/
        default:
            return "no label";
    }
}

/// A functor to access the surfaces of a volume
template <typename scalar_t, detray::n_axis::label axis_label>
struct edge_getter {
    /// Call operator that forwards the neighborhood search call in a volume
    /// to a surface finder data structure
    template <typename group_t, typename index_t,
              typename... Args>
    DETRAY_HOST_DEVICE inline auto operator()([[maybe_unused]] const group_t &group,
                                              [[maybe_unused]] const index_t index) const {
        using accel_t = typename group_t::value_type;
        if constexpr (detray::detail::is_grid_v<accel_t>) {
            std::cout << "enter \n";
            const auto grid = group[index];
            const auto axis_bin_edges = grid.axes().template get_axis<axis_label>().bin_edges();
            std::vector<scalar_t> edges;
            std::copy(axis_bin_edges.cbegin(), axis_bin_edges.cend(), std::back_inserter(edges));
            return edges;
        }
        std::cout << "no enter \n";
        return std::vector<scalar_t>{};
    }
};

template <typename detray::n_axis::label axis_label, typename detector_t, typename link_t>
auto bin_edges(const detector_t& detector, const link_t& link){
    using d_scalar_t = typename detector_t::scalar_type;
    return detector.surface_store().template visit<edge_getter<d_scalar_t, axis_label>>(link);
}

template <typename accel_ids_t, typename link_t>
actsvg::proto::grid::type get_grid_type(const link_t& link)
{
    std::cout << std::to_string(static_cast<int>(link.id())) + "\n";
    switch (link.id()){
            case accel_ids_t::e_cylinder2_grid:
                std::cout << "e_cylinder2_grid \n";
                return actsvg::proto::grid::e_x_y;
            case accel_ids_t::e_disc_grid:
                std::cout << "e_r_phi \n";
                return actsvg::proto::grid::e_r_phi;
            case accel_ids_t::e_brute_force:
                std::cout << "e_brute_force \n";
                return actsvg::proto::grid::e_x_y;
        }
    std::cout << "??? \n";
    return actsvg::proto::grid::e_x_y;
}

// Calculating r under the assumption that for edges_rphi then (bin_edge_min - bin_edge_max) / (2*pi).
template <typename d_scalar_t>
auto r_phi_split(const std::vector<d_scalar_t>& edges_rphi){
    d_scalar_t pi2 = (2.f * detray::constant<d_scalar_t>().pi);
    const auto [min, max] = std::minmax_element(edges_rphi.cbegin(), edges_rphi.cend());
    auto r = (*max - *min) / pi2;
    std::vector<d_scalar_t> edges_phi;
    std::transform(edges_rphi.cbegin(), edges_rphi.cend(), std::back_inserter(edges_phi), [r](d_scalar_t e){return e/r;});
    return std::tuple(edges_phi, r);
}

template <typename detector_t, typename link_t, typename view_t>
auto cylinder2_grid_type_and_edges(const detector_t& detector, const link_t& link, const view_t&){

    auto edges_rphi = bin_edges<detray::n_axis::label::e_rphi>(detector, link);
    auto edges_z = bin_edges<detray::n_axis::label::e_cyl_z>(detector, link);
    /*auto [edges_phi, r] = r_phi_split(edges_rphi);
    std::vector edges_r{0.f, r}

    if (std::is_same_v<view_t, actsvg::views::x_y>()){
        return std::tuple(actsvg::proto::grid::e_r_phi, edges_r, edges_phi);
    }
    if (std::is_same_v<view_t, actsvg::views::z_r>()){
        return std::tuple(actsvg::proto::grid::e_x_y, edges_z, edges_r);
    }
    if (std::is_same_v<view_t, actsvg::views::z_phi>()){
        return std::tuple(actsvg::proto::grid::e_x_y, edges_z, edges_phi);
    }*/
    if (std::is_same_v<view_t, typename actsvg::views::z_rphi>){
        return std::tuple(actsvg::proto::grid::e_x_y, edges_z, edges_rphi);
    }
    //throw std::domain_error("View must x_y, z_r, z_phi, or z_rphi");
    using scalar_t = typename detector_t::scalar_type;
    return std::tuple(actsvg::proto::grid::e_x_y, std::vector<scalar_t>{}, std::vector<scalar_t>{});
}

template <typename detector_t, typename link_t, typename view_t>
auto disc_grid_type_and_edges(const detector_t& detector, const link_t& link, const view_t&){

    auto edges_r = bin_edges<detray::n_axis::label::e_r>(detector, link);
    auto edges_phi = bin_edges<detray::n_axis::label::e_phi>(detector, link);

    if (std::is_same_v<view_t, typename actsvg::views::x_y>){
        return std::tuple(actsvg::proto::grid::e_r_phi, edges_r, edges_phi);
    }
    /*if (std::is_same_v<view_t, actsvg::views::z_r>()){
        return std::tuple(actsvg::proto::grid::e_x_y, edges_z, edges_r);
    }
    if (std::is_same_v<view_t, actsvg::views::z_phi>()){
        return std::tuple(actsvg::proto::grid::e_x_y, edges_z, edges_phi);
    }
    if (std::is_same_v<view_t, actsvg::views::z_rphi>()){
        return std::tuple(actsvg::proto::grid::e_x_y, edges_z, edges_rphi);
    }*/
    //throw std::domain_error("View must x_y, z_r, z_phi, or z_rphi");
    using scalar_t = typename detector_t::scalar_type;
    return std::tuple(actsvg::proto::grid::e_x_y, std::vector<scalar_t>{}, std::vector<scalar_t>{});
}

template <typename accel_ids_t, typename detector_t, typename link_t, typename view_t>
auto get_type_and_axes(const detector_t& detector, const link_t& link, const view_t& view)
{
    switch (link.id()){
        case accel_ids_t::e_cylinder2_grid: {
            return cylinder2_grid_type_and_edges(detector, link, view);
        }
        case accel_ids_t::e_disc_grid: {
            return disc_grid_type_and_edges(detector, link, view);
        }
        default:
        {
            using scalar_t = typename detector_t::scalar_type;
            return std::tuple(actsvg::proto::grid::e_x_y, std::vector<scalar_t>{}, std::vector<scalar_t>{});
        }
    }
}

template <typename a_scalar_t, typename detector_t, typename view_t>
auto grid(const detector_t& detector, const std::size_t index, const view_t& view)
{
    using d_scalar_t = typename detector_t::scalar_type;
    using geo_object_ids = typename detector_t::geo_obj_ids;
    using accel_ids = typename detector_t::sf_finders::id;

    const auto vol_desc = detector.volumes()[index];
    const auto link = vol_desc.template link<geo_object_ids::e_sensitive>();
    actsvg::proto::grid p_grid;

    if (not link.is_invalid()) {
        const auto [type, edges0, edges1] = get_type_and_axes<accel_ids>(detector, link, view);
        p_grid._type = type;
        std::transform(edges0.cbegin(), edges0.cend(), std::back_inserter(p_grid._edges_0), [](d_scalar_t v){ return static_cast<a_scalar_t>(v); });
        std::transform(edges1.cbegin(), edges1.cend(), std::back_inserter(p_grid._edges_1), [](d_scalar_t v){ return static_cast<a_scalar_t>(v); });

        for (auto i: p_grid._edges_0)
            std::cout << std::to_string(i) << ' ';
        std::cout << "\n NEXT: \n";
        for (auto i: p_grid._edges_1)
            std::cout << std::to_string(i) << ' ';
        std::cout << "\n";
    }
    return p_grid;
}

template <typename detector_t, typename view_t>
auto draw_grid(const std::string& identification, const detector_t& detector, const std::size_t index, const view_t& view){
    auto p_grid = grid<actsvg::scalar>(detector, index, view);
    p_grid._stroke._width = 5.f;
    return actsvg::display::grid(identification, p_grid);
}

int main(int, char**) {
    
    // Creating the detector and geomentry context.
    using detector_t = detray::detector<detray::toy_metadata<>>;
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::create_toy_geometry(host_mr, 4, 3);
    detector_t::geometry_context context{};

    // Creating the view.
    const actsvg::views::z_r view;

    // Creating the svg generator for the detector.
    detray::svgtools::illustrator il{det, names, true};

    std::cout << "volumes size: " + std::to_string(det.volumes().size()) + "\n";

    std::vector<std::size_t> indices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                        10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    
    for (const auto i : indices){
        
        std::string name = "volume" + std::to_string(i) + "_grid";
        std::cout << "checking " + name + "\n";
        const auto grid_svg = draw_grid(name, det, i, view);
        const auto volume_svg = il.draw_volume("volume", context, i, view);
        
        detray::svgtools::write_svg(name + ".svg", {volume_svg, grid_svg});
        std::cout << "\n";
    }
}