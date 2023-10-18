/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/grids/axis.hpp"

// Actsvg include(s)
#include "actsvg/proto/grid.hpp"

// System include(s)
#include <algorithm>
#include <optional>
#include <string>
#include <vector>

namespace detray::svgtools::conversion {

/// A functor to access the bin edges.
template <typename scalar_t, detray::n_axis::label axis_label>
struct edge_getter {
    template <typename group_t, typename index_t, typename... Args>
    DETRAY_HOST_DEVICE inline auto operator()(
        [[maybe_unused]] const group_t& group,
        [[maybe_unused]] const index_t index) const {
        using accel_t = typename group_t::value_type;
        if constexpr (detray::detail::is_grid_v<accel_t>) {
            const auto grid = group[index];
            const auto axis_bin_edges =
                grid.axes().template get_axis<axis_label>().bin_edges();
            std::vector<scalar_t> edges;
            std::copy(axis_bin_edges.cbegin(), axis_bin_edges.cend(),
                      std::back_inserter(edges));
            return edges;
        }
        return std::vector<scalar_t>{};
    }
};

/// @return the bin edges of the grid.
template <typename detray::n_axis::label axis_label, typename detector_t,
          typename link_t>
auto bin_edges(const detector_t& detector, const link_t& link) {
    using d_scalar_t = typename detector_t::scalar_type;
    return detector.surface_store()
        .template visit<edge_getter<d_scalar_t, axis_label>>(link);
}

// Calculating r under the assumption that the cylinder grid is closed in phi:
// (bin_edge_min - bin_edge_max) / (2*pi).
template <typename d_scalar_t>
auto r_phi_split(const std::vector<d_scalar_t>& edges_rphi) {
    const auto r = edges_rphi.back() * detray::constant<d_scalar_t>::inv_pi;
    std::vector<d_scalar_t> edges_phi;
    std::transform(edges_rphi.cbegin(), edges_rphi.cend(),
                   std::back_inserter(edges_phi),
                   [r](d_scalar_t e) { return e / r; });
    return std::tuple(edges_phi, r);
}

/// @returns the actsvg grid type and edge values for a detray cylinder
/// grid.
template <typename detector_t, typename link_t, typename view_t>
auto cylinder2_grid_type_and_edges(const detector_t& detector,
                                   const link_t& link, const view_t&) {
    assert(link.id() == detector_t::sf_finders::id::e_cylinder2_grid);
    auto edges_rphi = bin_edges<detray::n_axis::label::e_rphi>(detector, link);
    auto edges_z = bin_edges<detray::n_axis::label::e_cyl_z>(detector, link);
    auto [edges_phi, r] = r_phi_split(edges_rphi);
    std::vector edges_r{r, r};

    if (std::is_same_v<view_t, actsvg::views::x_y>) {
        return std::tuple(actsvg::proto::grid::e_r_phi, edges_r, edges_phi);
    }
    if (std::is_same_v<view_t, actsvg::views::z_r>) {
        return std::tuple(actsvg::proto::grid::e_x_y, edges_z, edges_r);
    }
    if (std::is_same_v<view_t, actsvg::views::z_phi>) {
        return std::tuple(actsvg::proto::grid::e_x_y, edges_z, edges_phi);
    }
    if (std::is_same_v<view_t, typename actsvg::views::z_rphi>) {
        return std::tuple(actsvg::proto::grid::e_x_y, edges_z, edges_rphi);
    }
    using scalar_t = typename detector_t::scalar_type;
    return std::tuple(actsvg::proto::grid::e_x_y, std::vector<scalar_t>{},
                      std::vector<scalar_t>{});
}

/// @returns the actsvg grid type and edge values for a detray disc grid.
template <typename detector_t, typename link_t, typename view_t>
auto disc_grid_type_and_edges(const detector_t& detector, const link_t& link,
                              const view_t&) {
    assert(link.id() == detector_t::sf_finders::id::e_disc_grid);
    auto edges_r = bin_edges<detray::n_axis::label::e_r>(detector, link);
    auto edges_phi = bin_edges<detray::n_axis::label::e_phi>(detector, link);

    if (std::is_same_v<view_t, typename actsvg::views::x_y>) {
        return std::tuple(actsvg::proto::grid::e_r_phi, edges_r, edges_phi);
    }
    using scalar_t = typename detector_t::scalar_type;
    return std::tuple(actsvg::proto::grid::e_x_y, std::vector<scalar_t>{},
                      std::vector<scalar_t>{});
}

/// @returns the detray grids respective actsvg grid type and edge
/// values.
template <typename detector_t, typename link_t,
          typename view_t>
auto get_type_and_axes(const detector_t& detector, const link_t& link,
                       const view_t& view) {
    using accel_ids_t = typename detector_t::sf_finders::id;
    switch (link.id()) {
        case accel_ids_t::e_cylinder2_grid: {
            return cylinder2_grid_type_and_edges(detector, link, view);
        }
        case accel_ids_t::e_disc_grid: {
            return disc_grid_type_and_edges(detector, link, view);
        }
        default: {
            using scalar_t = typename detector_t::scalar_type;
            return std::tuple(actsvg::proto::grid::e_x_y,
                              std::vector<scalar_t>{}, std::vector<scalar_t>{});
        }
    }
}

/// @brief Converts a detray grid to a actsvg proto grid.
/// @param detector the detector
/// @param index the index of the grid's volume
/// @param view the view
/// @returns a proto grid
template <typename a_scalar_t, typename detector_t, typename view_t>
std::optional<actsvg::proto::grid> grid(const detector_t& detector,
                                        const std::size_t index,
                                        const view_t& view) {
    using d_scalar_t = typename detector_t::scalar_type;
    using geo_object_ids = typename detector_t::geo_obj_ids;
    using accel_ids = typename detector_t::sf_finders::id;

    const auto vol_desc = detector.volumes()[index];
    const auto link = vol_desc.template link<geo_object_ids::e_sensitive>();
    actsvg::proto::grid p_grid;

    if (not link.is_invalid()) {
        const auto [type, edges0, edges1] =
            get_type_and_axes(detector, link, view);
        p_grid._type = type;
        std::transform(edges0.cbegin(), edges0.cend(),
                       std::back_inserter(p_grid._edges_0),
                       [](d_scalar_t v) { return static_cast<a_scalar_t>(v); });
        std::transform(edges1.cbegin(), edges1.cend(),
                       std::back_inserter(p_grid._edges_1),
                       [](d_scalar_t v) { return static_cast<a_scalar_t>(v); });
        return p_grid;
    }
    return {};
}

}  // namespace detray::svgtools::conversion