/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/grid_axis.hpp"
#include "detray/definitions/units.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"
#include "detray/plugins/svgtools/utils/surface_kernels.hpp"

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
    using axis_label = detray::axis::label;

    auto edges_phi = grid.template get_axis<axis_label::e_rphi>().bin_edges();
    auto edges_z = grid.template get_axis<axis_label::e_cyl_z>().bin_edges();
    // Unknown for 2D cylinder
    dvector<scalar_t> edges_r{};

    if constexpr (std::is_same_v<view_t, actsvg::views::x_y>) {
        return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_r_phi,
                          edges_r, edges_phi);
    }
    if constexpr (std::is_same_v<view_t, actsvg::views::z_r>) {
        return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_x_y,
                          edges_z, edges_r);
    }
    if constexpr (std::is_same_v<view_t, actsvg::views::z_phi> ||
                  std::is_same_v<view_t, typename actsvg::views::z_rphi>) {
        return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_z_phi,
                          edges_z, edges_phi);
    }

    return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_x_y,
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
    using axis_label = detray::axis::label;

    auto edges_r = grid.template get_axis<axis_label::e_r>().bin_edges();
    auto edges_phi = grid.template get_axis<axis_label::e_phi>().bin_edges();

    if constexpr (std::is_same_v<view_t, typename actsvg::views::x_y>) {
        return std::tuple(grid_type::e_endcap, actsvg::proto::grid::e_r_phi,
                          edges_r, edges_phi);
    }

    return std::tuple(grid_type::e_endcap, actsvg::proto::grid::e_x_y, edges_r,
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

        return std::tuple(grid_type::e_unknown, actsvg::proto::grid::e_x_y,
                          dvector<scalar_t>{}, dvector<scalar_t>{});
    }
};

/// A functor to access the bins and get the associated surface indices
struct bin_association_getter {

    template <typename group_t, typename index_t, typename volume_t>
    DETRAY_HOST_DEVICE std::vector<std::vector<std::size_t>> operator()(
        [[maybe_unused]] const group_t& group,
        [[maybe_unused]] const index_t index,
        [[maybe_unused]] const volume_t& vol_desc,
        [[maybe_unused]] const std::array<dindex, 2>& search_window) const {

        using accel_t = typename group_t::value_type;

        if constexpr (detray::detail::is_grid_v<accel_t>) {

            using transform3_t =
                typename accel_t::local_frame_type::transform3_type;
            using scalar_t = typename transform3_t::scalar_type;
            using point2_t = typename transform3_t::point2;

            // The sheet display only works for 2-dimensional grids
            if constexpr (accel_t::dim != 2u) {
                std::cout
                    << "WARNIGN: Only 2D grids can be displayed as actvg sheets"
                    << std::endl;
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
                        // actsvg expects the sensitive surfaces to be numbered
                        // starting from zero
                        dindex offset{vol_desc.template sf_link<
                            surface_id::e_sensitive>()[0]};
                        entries.push_back(sf_desc.index() - offset);
                    }

                    // Remove duplicates
                    std::sort(entries.begin(), entries.end());
                    auto last = std::unique(entries.begin(), entries.end());
                    entries.erase(last, entries.end());

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

    const auto vol_desc = detector.volume(index);
    const auto link =
        vol_desc.template accel_link<geo_object_ids::e_sensitive>();
    actsvg::proto::grid p_grid;

    if (not link.is_invalid()) {
        auto [gr_type, view_type, edges0, edges1] =
            detector.accelerator_store()
                .template visit<detail::type_and_edge_getter<scalar_t>>(link,
                                                                        view);
        p_grid._type = view_type;

        // Find the correct grid radius
        if (gr_type == detail::grid_type::e_barrel) {
            // Get the the radii of the volume portal surfaces
            std::vector<scalar_t> radii{};
            const auto vol = detector_volume{detector, vol_desc};

            // Passive surfaces could be in the brute force finder, but no
            // sensitive surfaces, since the volume has a grid. Their radii are,
            // however, always within the interval of the portal radii
            for (const auto& pt_desc : vol.portals()) {
                auto r = detector.mask_store()
                             .template visit<
                                 detray::svgtools::utils::outer_radius_getter>(
                                 pt_desc.mask());
                if (r.has_value()) {
                    radii.push_back(*r);
                }
            }

            scalar_t inner_r = *std::min_element(radii.begin(), radii.end());
            scalar_t outer_r = *std::max_element(radii.begin(), radii.end());

            p_grid._reference_r =
                0.5f * static_cast<actsvg::scalar>(inner_r + outer_r);

            // Add the cylinder radius to the axis binning
            if constexpr (std::is_same_v<view_t, actsvg::views::x_y>) {
                if (edges0.empty()) {
                    edges0 = {p_grid._reference_r, p_grid._reference_r};
                }
            }
            if constexpr (std::is_same_v<view_t, actsvg::views::z_r>) {
                if (edges1.empty()) {
                    edges1 = {p_grid._reference_r, p_grid._reference_r};
                }
            }

        } else if (gr_type == detail::grid_type::e_endcap) {
            // An axis is always sorted
            p_grid._reference_r = static_cast<actsvg::scalar>(edges0.back());
        }

        // Transform cylinder grid to rphi edges, if rphi view is requested
        if constexpr (std::is_same_v<view_t, typename actsvg::views::z_rphi>) {
            if (gr_type == detail::grid_type::e_barrel) {
                for (auto& e : edges1) {
                    e *= p_grid._reference_r;
                }
            }
        }

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
    const std::array<dindex, 2>& search_window = {2u, 2u}) {

    using geo_object_ids = typename detector_t::geo_obj_ids;

    const auto& vol_desc = det.volume(vol.index());
    const auto& link =
        vol_desc.template accel_link<geo_object_ids::e_sensitive>();

    if (not link.is_invalid()) {
        return det.accelerator_store()
            .template visit<detail::bin_association_getter>(link, vol_desc,
                                                            search_window);
    }

    return {};
}

}  // namespace detray::svgtools::conversion
