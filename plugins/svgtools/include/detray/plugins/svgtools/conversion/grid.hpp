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
#include "detray/plugins/svgtools/styling/styling.hpp"

// Actsvg include(s)
#include "actsvg/proto/grid.hpp"

// System include(s)
#include <algorithm>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace detray::svgtools::conversion {

namespace detail {

enum grid_type : std::uint8_t { e_barrel = 0, e_endcap = 1, e_unknown = 2 };

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
        return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_r_phi, r,
                          edges_r, edges_phi);
    }
    if constexpr (std::is_same_v<view_t, actsvg::views::z_r>) {
        return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_x_y, r,
                          edges_z, edges_r);
    }
    if constexpr (std::is_same_v<view_t, actsvg::views::z_phi>) {
        return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_z_phi, r,
                          edges_z, edges_phi);
    }
    if constexpr (std::is_same_v<view_t, typename actsvg::views::z_rphi>) {
        return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_z_phi, r,
                          edges_z, edges_rphi);
    }

    return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_x_y, r,
                      dvector<scalar_t>{}, dvector<scalar_t>{});
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
        return std::tuple(grid_type::e_endcap, actsvg::proto::grid::e_r_phi, r,
                          edges_r, edges_phi);
    }

    return std::tuple(grid_type::e_endcap, actsvg::proto::grid::e_x_y, r,
                      dvector<scalar_t>{}, dvector<scalar_t>{});
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

        return std::tuple(grid_type::e_unknown, actsvg::proto::grid::e_x_y,
                          detray::detail::invalid_value<scalar_t>(),
                          dvector<scalar_t>{}, dvector<scalar_t>{});
    }
};

/// A functor to access the bins and get the associated surface indices
struct bin_association_getter {

    template <typename group_t, typename index_t>
    DETRAY_HOST_DEVICE std::vector<std::vector<std::size_t>> operator()(
        [[maybe_unused]] const group_t& group,
        [[maybe_unused]] const index_t index, [[maybe_unused]] dindex offset,
        [[maybe_unused]] const std::array<dindex, 2>& search_window) const {

        using accel_t = typename group_t::value_type;

        if constexpr (detray::detail::is_grid_v<accel_t>) {

            using transform3_t =
                typename accel_t::local_frame_type::transform3_type;
            using scalar_t = typename transform3_t::scalar_type;
            using point2_t = typename transform3_t::point2;

            // The sheet display only works for 2-dimensional grids
            if constexpr (accel_t::Dim != 2u) {
                return {};
            }

            const accel_t grid = group[index];
            const std::size_t n_bins{grid.nbins()};

            std::vector<std::vector<std::size_t>> bin_assoc;
            bin_assoc.reserve(n_bins);

            // Create the bin associations
            auto edges0 = grid.template get_axis<0>().bin_edges();
            auto edges1 = grid.template get_axis<1>().bin_edges();

            // In the svg convention the phi axis has to be the second axis to
            // loop over
            constexpr bool is_cyl{
                std::is_same_v<typename accel_t::local_frame_type,
                               detray::cylindrical2<transform3_t>>};
            if constexpr (is_cyl) {
                edges0.swap(edges1);
            }

            for (std::size_t i = 1u; i < edges0.size(); ++i) {
                scalar_t p0 = 0.5f * (edges0[i] + edges0[i - 1]);

                for (std::size_t j = 1u; j < edges1.size(); ++j) {
                    scalar_t p1 = 0.5f * (edges1[j] + edges1[j - 1]);

                    // Create the bin center position estimates
                    point2_t bin_center{p0, p1};
                    if constexpr (is_cyl) {
                        bin_center = {p1, p0};
                    }

                    // Get all the bin entries and calculate the loc index
                    std::vector<std::size_t> entries;

                    for (const auto& sf_desc :
                         grid.search(bin_center, search_window)) {
                        entries.push_back(sf_desc.index() - offset);
                    }

                    bin_assoc.push_back(std::move(entries));
                }
            }

            return bin_assoc;
        }

        return {};
    }
};

}  // namespace detail

/// @brief Converts a detray grid to a actsvg proto grid.
///
/// @param detector the detector
/// @param index the index of the grid's volume
/// @param view the view
/// @param style the style settings
///
/// @returns a proto grid
template <typename detector_t, typename view_t>
auto grid(const detector_t& detector, const dindex index, const view_t& view,
          const styling::grid_style& style =
              styling::tableau_colorblind::grid_style) {

    using scalar_t = typename detector_t::scalar_type;
    using geo_object_ids = typename detector_t::geo_obj_ids;

    const auto& vol_desc = detector.volumes()[index];
    const auto& link = vol_desc.template link<geo_object_ids::e_sensitive>();
    actsvg::proto::grid p_grid;

    if (not link.is_invalid()) {
        const auto [gr_type, view_type, r, edges0, edges1] =
            detector.accelerator_store()
                .template visit<detail::type_and_edge_getter<scalar_t>>(link,
                                                                        view);

        p_grid._type = view_type;
        p_grid._reference_r = static_cast<actsvg::scalar>(r);

        std::transform(
            edges0.cbegin(), edges0.cend(), std::back_inserter(p_grid._edges_0),
            [](scalar_t v) { return static_cast<actsvg::scalar>(v); });
        std::transform(
            edges1.cbegin(), edges1.cend(), std::back_inserter(p_grid._edges_1),
            [](scalar_t v) { return static_cast<actsvg::scalar>(v); });

        svgtools::styling::apply_style(p_grid, style);

        return std::tuple(std::optional<actsvg::proto::grid>{p_grid}, gr_type);
    }

    return std::tuple(std::optional<actsvg::proto::grid>{},
                      detail::grid_type::e_unknown);
}

/// @brief Get the surfaces indices that are registered in the bin neighborhoods
///
/// @param detector the detector
/// @param vol the volume that holds the grid and surfaces
/// @param offset transform a global surface index to a local one for the volume
///
/// @returns a vector of surface indices per neighborhood
template <typename detector_t>
std::vector<std::vector<std::size_t>> get_bin_association(
    const detector_t& det, const detray::detector_volume<detector_t>& vol,
    std::size_t offset = 0u,
    const std::array<dindex, 2>& search_window = {2u, 2u}) {

    using geo_object_ids = typename detector_t::geo_obj_ids;

    const auto& vol_desc = det.volumes()[vol.index()];
    const auto& link = vol_desc.template link<geo_object_ids::e_sensitive>();

    if (not link.is_invalid()) {
        return det.accelerator_store()
            .template visit<detail::bin_association_getter>(
                link, static_cast<dindex>(offset), search_window);
    }

    return {};
}

}  // namespace detray::svgtools::conversion
