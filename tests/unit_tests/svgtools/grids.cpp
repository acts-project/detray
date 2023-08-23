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
template <typename scalar_t, detray::n_axis::label axis1, detray::n_axis::label axis2>
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
            const auto axis0_bin_edges = grid.axes().template get_axis<axis1>().bin_edges();
            const auto axis1_bin_edges = grid.axes().template get_axis<axis2>().bin_edges();
            std::vector<scalar_t> edges0;
            std::vector<scalar_t> edges1;
            std::copy(axis0_bin_edges.cbegin(), axis0_bin_edges.cend(), std::back_inserter(edges0));
            std::copy(axis1_bin_edges.cbegin(), axis1_bin_edges.cend(), std::back_inserter(edges1));
            return std::tuple(edges0, edges1);
        }
        std::cout << "no enter \n";
        return std::tuple(std::vector<scalar_t>{}, std::vector<scalar_t>{});
    }
};
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

template <typename accel_ids_t, typename detector_t, typename link_t>
auto get_type_and_axes(const detector_t& detector, const link_t& link, const actsvg::views::x_y&)
{
    using d_scalar_t = typename detector_t::scalar_type;
    switch (link.id()){
        case accel_ids_t::e_cylinder2_grid: {
            auto [edges_rphi, edges_z] = detector.surface_store().template visit<edge_getter<d_scalar_t, detray::n_axis::label::e_rphi, detray::n_axis::label::e_cyl_z>>(link);
            d_scalar_t pi2 = (1.f * detray::constant<d_scalar_t>().pi); //Double check this...
            auto r = *std::max_element(edges_rphi.cbegin(), edges_rphi.cend()) / pi2;
            std::vector<d_scalar_t> edges_phi;
            std::transform(edges_rphi.cbegin(), edges_rphi.cend(), std::back_inserter(edges_phi), [r](d_scalar_t e){return e/r;});
            return std::tuple(actsvg::proto::grid::e_r_phi, std::vector{0.f, r}, edges_phi);
        }
        default:{
            // Brute forcer and disc grid:
            // (Brute forcer does not display anything)
            auto [edges_r, edges_phi] = detector.surface_store().template visit<edge_getter<d_scalar_t, detray::n_axis::label::e_r, detray::n_axis::label::e_phi>>(link);
            return std::tuple(actsvg::proto::grid::e_r_phi, edges_r, edges_phi);
        }
    }
}

template <typename a_scalar_t, typename detector_t>
auto grid(const detector_t& detector, const std::size_t index)
{
    using d_scalar_t = typename detector_t::scalar_type;
    using geo_object_ids = typename detector_t::geo_obj_ids;
    using accel_ids = typename detector_t::sf_finders::id;

    const auto vol_desc = detector.volumes()[index];
    const auto link = vol_desc.template link<geo_object_ids::e_sensitive>();

    actsvg::proto::grid p_grid;

    if (not link.is_invalid()) {
        const auto [type, edges0, edges1] = get_type_and_axes<accel_ids>(detector, link, actsvg::views::x_y{});
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

template <typename detector_t>
auto draw_grid(const std::string& identification, const detector_t& detector, const std::size_t index){
    auto p_grid = grid<actsvg::scalar>(detector, index);
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
    const actsvg::views::x_y view;

    // Creating the svg generator for the detector.
    detray::svgtools::illustrator il{det, names, true};

    std::cout << "volumes size: " + std::to_string(det.volumes().size()) + "\n";

    std::vector<std::size_t> indices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                        10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    
    for (const auto i : indices){
        
        std::string name = "volume" + std::to_string(i) + "_grid";
        std::cout << "checking " + name + "\n";
        const auto grid_svg = draw_grid(name, det, i);
        const auto volume_svg = il.draw_volume("volume", context, i, view);
        
        detray::svgtools::write_svg(name + ".svg", {volume_svg, grid_svg});
        std::cout << "\n";
    }
}