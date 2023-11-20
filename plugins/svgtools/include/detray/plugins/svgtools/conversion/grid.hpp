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

namespace detail {

/// @returns the phi edges and radius of an r-phi axis
template <typename scalar_t>
inline auto r_phi_split(const dvector<scalar_t>& edges_rphi) {

    // Calculating r under the assumption that the cylinder grid is closed in
    // phi: (bin_edge_min - bin_edge_max) / (2*pi).
    const auto r = edges_rphi.back() * detray::constant<scalar_t>::inv_pi;
    dvector<scalar_t> edges_phi;

    std::transform(edges_rphi.cbegin(), edges_rphi.cend(),
                   std::back_inserter(edges_phi),
                   [r](scalar_t e) { return e / r; });

    return std::tuple(edges_phi, r);
}

/// @returns the actsvg grid type and edge values for a detray 2D cylinder grid.
template <
    typename grid_t, typename view_t,
    std::enable_if_t<
        std::is_same_v<typename grid_t::local_frame_type,
                       detray::cylindrical2<
                           typename grid_t::local_frame_type::transform3_type>>,
        bool> = true>
inline auto grid_type_and_edges(const grid_t& grid, const view_t&) {

    using scalar_t = typename grid_t::local_frame_type::scalar_type;
    using axis_label = detray::n_axis::label;

    auto edges_rphi = grid.template get_axis<axis_label::e_rphi>().bin_edges();
    auto edges_z = grid.template get_axis<axis_label::e_cyl_z>().bin_edges();
    auto [edges_phi, r] = r_phi_split(edges_rphi);
    dvector<scalar_t> edges_r{r, r};

    if constexpr (std::is_same_v<view_t, actsvg::views::x_y>) {
        return std::tuple(actsvg::proto::grid::e_r_phi, r, edges_r, edges_phi);
    }
    if constexpr (std::is_same_v<view_t, actsvg::views::z_r>) {
        return std::tuple(actsvg::proto::grid::e_x_y, r, edges_z, edges_r);
    }
    if constexpr (std::is_same_v<view_t, actsvg::views::z_phi>) {
        return std::tuple(actsvg::proto::grid::e_z_phi, r, edges_z, edges_phi);
    }
    if constexpr (std::is_same_v<view_t, typename actsvg::views::z_rphi>) {
        return std::tuple(actsvg::proto::grid::e_z_phi, r, edges_z, edges_rphi);
    }

    return std::tuple(actsvg::proto::grid::e_x_y, r, dvector<scalar_t>{},
                      dvector<scalar_t>{});
}

/// @returns the actsvg grid type and edge values for a detray disc grid.
template <
    typename grid_t, typename view_t,
    std::enable_if_t<
        std::is_same_v<
            typename grid_t::local_frame_type,
            detray::polar2<typename grid_t::local_frame_type::transform3_type>>,
        bool> = true>
inline auto grid_type_and_edges(const grid_t& grid, const view_t&) {

    using scalar_t = typename grid_t::local_frame_type::scalar_type;
    using axis_label = detray::n_axis::label;

    auto edges_r = grid.template get_axis<axis_label::e_r>().bin_edges();
    auto edges_phi = grid.template get_axis<axis_label::e_phi>().bin_edges();
    // The axes are always sorted
    scalar_t r{edges_r.back()};

    if constexpr (std::is_same_v<view_t, typename actsvg::views::x_y>) {
        return std::tuple(actsvg::proto::grid::e_r_phi, r, edges_r, edges_phi);
    }

    return std::tuple(actsvg::proto::grid::e_x_y, r, dvector<scalar_t>{},
                      dvector<scalar_t>{});
}

/// A functor to access the type and bin edges of a grid.
template <typename scalar_t>
struct type_and_edge_getter {

    template <typename group_t, typename index_t, typename view_t>
    DETRAY_HOST_DEVICE inline auto operator()(
        [[maybe_unused]] const group_t& group,
        [[maybe_unused]] const index_t index,
        [[maybe_unused]] const view_t& view) const {

        using accel_t = typename group_t::value_type;

        if constexpr (detray::detail::is_grid_v<accel_t>) {
            return grid_type_and_edges(group[index], view);
        }

        return std::tuple(actsvg::proto::grid::e_x_y,
                          detray::detail::invalid_value<scalar_t>(),
                          dvector<scalar_t>{}, dvector<scalar_t>{});
    }
};

}  // namespace detail

/// @brief Converts a detray grid to a actsvg proto grid.
///
/// @param detector the detector
/// @param index the index of the grid's volume
/// @param view the view
///
/// @returns a proto grid
template <typename a_scalar_t, typename detector_t, typename view_t>
std::optional<actsvg::proto::grid> grid(const detector_t& detector,
                                        const std::size_t index,
                                        const view_t& view) {

    using d_scalar_t = typename detector_t::scalar_type;
    using geo_object_ids = typename detector_t::geo_obj_ids;

    const auto& vol_desc = detector.volumes()[index];
    const auto& link = vol_desc.template link<geo_object_ids::e_sensitive>();
    actsvg::proto::grid p_grid;

    if (not link.is_invalid()) {
        const auto [type, r, edges0, edges1] =
            detector.accelerator_store()
                .template visit<detail::type_and_edge_getter<d_scalar_t>>(link,
                                                                          view);

        p_grid._type = type;
        p_grid._reference_r = static_cast<a_scalar_t>(r);

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
