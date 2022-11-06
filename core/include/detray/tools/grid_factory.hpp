/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/masks/masks.hpp"
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/grid.hpp"
#include "detray/surface_finders/grid/grid_collection.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <cassert>
#include <iostream>
#include <string>

namespace detray {

namespace detail {

/// @brief Helper type to assemble an multi-axis from bounds tuple and binnings
template <bool is_owning, typename containers, typename local_frame, typename,
          typename>
struct multi_axis_assembler;

/// @brief Specialized struct to extract axis bounds from a tuple
template <bool is_owning, typename containers, typename local_frame,
          typename... axis_bounds, typename... binning_ts>
struct multi_axis_assembler<is_owning, containers, local_frame,
                            dtuple<axis_bounds...>, dtuple<binning_ts...>> {

    static_assert(sizeof...(axis_bounds) == sizeof...(binning_ts),
                  "Number of axis bounds for this mask and given binning types "
                  "don't match!");

    using type =
        n_axis::multi_axis<is_owning, local_frame,
                           n_axis::single_axis<axis_bounds, binning_ts>...>;
};

}  // namespace detail

/// Typedef for easier construction @c multi_axis from mask shapes
template <typename shape_t, bool is_owning = true,
          typename containers = host_container_types,
          typename algebra_t = __plugin::transform3<detray::scalar>>
using coordinate_axes = typename detail::multi_axis_assembler<
    is_owning, containers,
    typename shape_t::template coordinate_type<algebra_t>,
    typename shape_t::axes::types,
    typename shape_t::axes::template binning<
        containers, typename algebra_t::scalar_type>>::type;

/// @brief Provides functionality to instantiate grids and grid collections.
///
/// @tparam value_t type of entries in the grid bins
/// @tparam serialzier_t  type of the serializer to the storage represenations
/// @tparam populator_impl_t how to fill the bin content
/// @tparam algebra_t the matrix/vector/point types to use
/// @tparam container_t the container types to use
template <typename value_t, template <std::size_t> class serializer_t,
          typename populator_impl_t,
          typename algebra_t = __plugin::transform3<detray::scalar>>
class grid_factory {

    public:
    // All grids are owning since they are used to fill the data
    static constexpr bool is_owning = true;

    using value_type = value_t;
    template <typename grid_shape_t>
    using grid_type =
        grid<grid_shape_t, value_type, serializer_t, populator_impl_t>;

    using scalar_type = typename algebra_t::scalar_type;
    template <typename T>
    using vector_type = host_container_types::template vector_type<T>;
    using algebra_type = algebra_t;

    grid_factory() = default;

    /// Takes the resource of the detector to allocate memory correctly
    explicit grid_factory(vecmem::memory_resource &resource)
        : m_resource(&resource) {}

    /// Get new grid collection
    template <typename grid_t>
    constexpr auto new_collection() {
        return grid_collection<typename grid_t::template type<!is_owning>>(
            m_resource);
    }

    /// Print grid - up to three dimensions
    /// @note will likely become obsolete with the actsvg implementation
    template <typename grid_t>
    auto to_string(const grid_t &) const noexcept -> void;

    //
    // annulus 2D
    //
    template <
        typename r_bounds = n_axis::closed<n_axis::label::e_r>,
        typename phi_bounds = n_axis::circular<>,
        typename r_binning = n_axis::regular<host_container_types, scalar_type>,
        typename phi_binning =
            n_axis::regular<host_container_types, scalar_type>>
    auto new_grid(const mask<annulus2D<>> &grid_bounds,
                  const std::array<std::size_t, 2UL> n_bins,
                  const std::array<std::vector<scalar_type>, 2UL> &bin_edges = {
                      {}}) const {
        // Axes boundaries and local indices
        using boundary = annulus2D<>::boundaries;
        using axes = annulus2D<>::axes<>;
        constexpr auto e_r_axis = static_cast<std::size_t>(axes::axis_loc0);
        constexpr auto e_phi_axis = static_cast<std::size_t>(axes::axis_loc1);
        auto b_values = grid_bounds.values();

        return new_grid<annulus2D<>>(
            {b_values[boundary::e_min_r], b_values[boundary::e_max_r],
             b_values[boundary::e_average_phi] -
                 b_values[boundary::e_min_phi_rel],
             b_values[boundary::e_average_phi] +
                 b_values[boundary::e_max_phi_rel]},
            {n_bins[e_r_axis], n_bins[e_phi_axis]},
            {bin_edges[e_r_axis], bin_edges[e_phi_axis]},
            std::tuple<r_bounds, phi_bounds>{},
            std::tuple<r_binning, phi_binning>{});
    }

    //
    // cuboid 3D
    //
    template <
        typename x_bounds = n_axis::closed<n_axis::label::e_x>,
        typename y_bounds = n_axis::closed<n_axis::label::e_y>,
        typename z_bounds = n_axis::closed<n_axis::label::e_z>,
        typename x_binning = n_axis::regular<host_container_types, scalar_type>,
        typename y_binning = n_axis::regular<host_container_types, scalar_type>,
        typename z_binning = n_axis::regular<host_container_types, scalar_type>>
    auto new_grid(const mask<cuboid3D> &grid_bounds,
                  const std::array<std::size_t, 3UL> n_bins,
                  const std::array<std::vector<scalar_type>, 3UL> &bin_edges = {
                      {}}) const {
        // Axes boundaries and local indices
        using boundary = cuboid3D::boundaries;
        using axes = cuboid3D::axes<>;
        constexpr auto e_x_axis = static_cast<std::size_t>(axes::axis_loc0);
        constexpr auto e_y_axis = static_cast<std::size_t>(axes::axis_loc1);
        constexpr auto e_z_axis = static_cast<std::size_t>(axes::axis_loc2);

        auto b_values = grid_bounds.values();

        return new_grid<cuboid3D>(
            {-b_values[boundary::e_half_x], b_values[boundary::e_half_x],
             -b_values[boundary::e_half_y], b_values[boundary::e_half_y],
             -b_values[boundary::e_half_z], b_values[boundary::e_half_z]},
            {n_bins[e_x_axis], n_bins[e_y_axis], n_bins[e_z_axis]},
            {bin_edges[e_x_axis], bin_edges[e_y_axis], bin_edges[e_z_axis]},
            std::tuple<x_bounds, y_bounds, z_bounds>{},
            std::tuple<x_binning, y_binning, z_binning>{});
    }

    //
    // cylinder 2D
    //
    template <
        typename rphi_bounds = n_axis::circular<n_axis::label::e_rphi>,
        typename z_bounds = n_axis::closed<n_axis::label::e_cyl_z>,
        typename rphi_binning =
            n_axis::regular<host_container_types, scalar_type>,
        typename z_binning = n_axis::regular<host_container_types, scalar_type>>
    auto new_grid(const mask<cylinder2D<>> &grid_bounds,
                  const std::array<std::size_t, 2UL> n_bins,
                  const std::array<std::vector<scalar_type>, 2UL> &bin_edges = {
                      {}}) const {
        // Axes boundaries and local indices
        using boundary = cylinder2D<>::boundaries;
        using axes = cylinder2D<>::axes<>;
        constexpr auto e_rphi_axis = static_cast<std::size_t>(axes::axis_loc0);
        constexpr auto e_z_axis = static_cast<std::size_t>(axes::axis_loc1);

        auto b_values = grid_bounds.values();

        return new_grid<cylinder2D<>>(
            {-static_cast<scalar_type>(M_PI) * b_values[boundary::e_r],
             static_cast<scalar_type>(M_PI) * b_values[boundary::e_r],
             b_values[boundary::e_n_half_z], b_values[boundary::e_p_half_z]},
            {n_bins[e_rphi_axis], n_bins[e_z_axis]},
            {bin_edges[e_rphi_axis], bin_edges[e_z_axis]},
            std::tuple<rphi_bounds, z_bounds>{},
            std::tuple<rphi_binning, z_binning>{});
    }

    //
    // cylinder 3D
    //
    template <
        typename r_bounds = n_axis::closed<n_axis::label::e_r>,
        typename phi_bounds = n_axis::circular<>,
        typename z_bounds = n_axis::closed<n_axis::label::e_z>,
        typename r_binning = n_axis::regular<host_container_types, scalar_type>,
        typename phi_binning =
            n_axis::regular<host_container_types, scalar_type>,
        typename z_binning = n_axis::regular<host_container_types, scalar_type>>
    auto new_grid(const mask<cylinder3D> &grid_bounds,
                  const std::array<std::size_t, 3UL> n_bins,
                  const std::array<std::vector<scalar_type>, 3UL> &bin_edges = {
                      {}}) const {
        // Axes boundaries and local indices
        using boundary = cylinder3D::boundaries;
        using axes = cylinder3D::axes<>;
        constexpr auto e_r_axis = static_cast<std::size_t>(axes::axis_loc0);
        constexpr auto e_phi_axis = static_cast<std::size_t>(axes::axis_loc1);
        constexpr auto e_z_axis = static_cast<std::size_t>(axes::axis_loc2);
        auto b_values = grid_bounds.values();

        return new_grid<cylinder3D>(
            {0.f, b_values[boundary::e_r], -static_cast<scalar_type>(M_PI),
             static_cast<scalar_type>(M_PI), -b_values[boundary::e_n_half_z],
             b_values[boundary::e_p_half_z]},
            {n_bins[e_r_axis], n_bins[e_phi_axis], n_bins[e_z_axis]},
            {bin_edges[e_r_axis], bin_edges[e_phi_axis], bin_edges[e_z_axis]},
            std::tuple<r_bounds, phi_bounds, z_bounds>{},
            std::tuple<r_binning, phi_binning, z_binning>{});
    }

    //
    // polar 2D
    //
    template <
        typename r_bounds = n_axis::closed<n_axis::label::e_r>,
        typename phi_bounds = n_axis::circular<>,
        typename r_binning = n_axis::regular<host_container_types, scalar_type>,
        typename phi_binning =
            n_axis::regular<host_container_types, scalar_type>>
    auto new_grid(const mask<ring2D<>> &grid_bounds,
                  const std::array<std::size_t, 2UL> n_bins,
                  const std::array<std::vector<scalar_type>, 2UL> &bin_edges = {
                      {}}) const {
        // Axes boundaries and local indices
        using boundary = ring2D<>::boundaries;
        using axes = ring2D<>::axes<>;
        constexpr auto e_r_axis = static_cast<std::size_t>(axes::axis_loc0);
        constexpr auto e_phi_axis = static_cast<std::size_t>(axes::axis_loc1);

        auto b_values = grid_bounds.values();

        return new_grid<ring2D<>>(
            {b_values[boundary::e_inner_r], b_values[boundary::e_outer_r],
             -static_cast<scalar_type>(M_PI), static_cast<scalar_type>(M_PI)},
            {n_bins[e_r_axis], n_bins[e_phi_axis]},
            {bin_edges[e_r_axis], bin_edges[e_phi_axis]},
            std::tuple<r_bounds, phi_bounds>{},
            std::tuple<r_binning, phi_binning>{});
    }

    //
    // rectangle 2D
    //
    template <
        typename x_bounds = n_axis::closed<n_axis::label::e_x>,
        typename y_bounds = n_axis::closed<n_axis::label::e_y>,
        typename x_binning = n_axis::regular<host_container_types, scalar_type>,
        typename y_binning = n_axis::regular<host_container_types, scalar_type>>
    auto new_grid(const mask<rectangle2D<>> &grid_bounds,
                  const std::array<std::size_t, 2UL> n_bins,
                  const std::array<std::vector<scalar_type>, 2UL> &bin_edges = {
                      {}}) const {
        // Axes boundaries and local indices
        using boundary = rectangle2D<>::boundaries;
        using axes = rectangle2D<>::axes<>;
        constexpr auto e_x_axis = static_cast<std::size_t>(axes::axis_loc0);
        constexpr auto e_y_axis = static_cast<std::size_t>(axes::axis_loc1);

        auto b_values = grid_bounds.values();

        return new_grid<rectangle2D<>>(
            {-b_values[boundary::e_half_x], b_values[boundary::e_half_x],
             -b_values[boundary::e_half_y], b_values[boundary::e_half_y]},
            {n_bins[e_x_axis], n_bins[e_y_axis]},
            {bin_edges[e_x_axis], bin_edges[e_y_axis]},
            std::tuple<x_bounds, y_bounds>{},
            std::tuple<x_binning, y_binning>{});
    }

    //
    // trapezoid 2D
    //
    template <
        typename x_bounds = n_axis::closed<n_axis::label::e_x>,
        typename y_bounds = n_axis::closed<n_axis::label::e_y>,
        typename x_binning = n_axis::regular<host_container_types, scalar_type>,
        typename y_binning = n_axis::regular<host_container_types, scalar_type>>
    auto new_grid(const mask<trapezoid2D<>> &grid_bounds,
                  const std::array<std::size_t, 2UL> n_bins,
                  const std::array<std::vector<scalar_type>, 2UL> &bin_edges = {
                      {}}) const {
        // Axes boundaries and local indices
        using boundary = trapezoid2D<>::boundaries;
        using axes = trapezoid2D<>::axes<>;
        constexpr auto e_x_axis = static_cast<std::size_t>(axes::axis_loc0);
        constexpr auto e_y_axis = static_cast<std::size_t>(axes::axis_loc1);

        auto b_values = grid_bounds.values();

        return new_grid<trapezoid2D<>>(
            {-b_values[boundary::e_half_length_1],
             b_values[boundary::e_half_length_1],
             -b_values[boundary::e_half_length_2],
             b_values[boundary::e_half_length_2]},
            {n_bins[e_x_axis], n_bins[e_y_axis]},
            {bin_edges[e_x_axis], bin_edges[e_y_axis]},
            std::tuple<x_bounds, y_bounds>{},
            std::tuple<x_binning, y_binning>{});
    }

    /// @brief Create and empty grid with fully initialized axes.
    ///
    /// @tparam grid_shape_t the shape of the resulting grid
    ///         (e.g. cylinder2D).
    /// @tparam e_bounds the bounds of the regular axes
    ///         (open vs. closed bounds).
    /// @tparam binnings the binning types of the axes
    ///         (regular vs. irregular)
    ///
    /// @param spans the span of the axis values for regular axes, otherwise
    ///              ignored.
    /// @param n_bins the number of bins for regular axes, otherwise ignored
    /// @param ax_bin_edges the explicit bin edges for irregular axes
    ///                     (lower bin edges + the the upper edge of the
    ///                     last bin), otherwise ignored.
    template <typename grid_shape_t, typename... bound_ts,
              typename... binning_ts>
    auto new_grid(const std::vector<scalar_type> spans,
                  const std::vector<std::size_t> n_bins,
                  const std::vector<std::vector<scalar_type>> &ax_bin_edges,
                  const std::tuple<bound_ts...>,
                  const std::tuple<binning_ts...>) const {

        static_assert(sizeof...(bound_ts) == sizeof...(binning_ts),
                      "number of axis bounds and binning types has to match");

        // Build the coordinate axes and the grid
        using axes_t = n_axis::multi_axis<
            is_owning,
            typename grid_shape_t::template local_frame_type<algebra_t>,
            n_axis::single_axis<bound_ts, binning_ts>...>;
        using bin_t = typename grid_type<axes_t>::bin_type;

        // Prepare data
        vector_type<dindex_range> axes_data{};
        vector_type<scalar_type> bin_edges{};

        // Call init for every axis
        unroll_axis_init<std::tuple<binning_ts...>>(
            spans, n_bins, ax_bin_edges, axes_data, bin_edges,
            std::make_index_sequence<axes_t::Dim>{});

        // Assemble the grid and return it
        axes_t axes(std::move(axes_data), std::move(bin_edges));

        vector_type<bin_t> bin_data{};
        bin_data.resize(axes.nbins()[0] * axes.nbins()[1],
                        populator_impl_t::template init<
                            typename grid_type<axes_t>::value_type>());

        return grid_type<axes_t>(std::move(bin_data), std::move(axes));
    }

    private:
    /// Initialize a single axis (either regular or irregular)
    /// @note change to template lambda as soon as it becomes available.
    template <std::size_t I, typename binnings>
    auto axis_init([[maybe_unused]] const std::vector<scalar_type> &spans,
                   [[maybe_unused]] const std::vector<std::size_t> &n_bins,
                   [[maybe_unused]] const std::vector<std::vector<scalar_type>>
                       &ax_bin_edges,
                   vector_type<dindex_range> &axes_data,
                   vector_type<scalar_type> &bin_edges) const {
        if constexpr (std::is_same_v<
                          std::tuple_element_t<I, binnings>,
                          n_axis::regular<host_container_types, scalar_type>>) {
            axes_data.push_back({bin_edges.size(), n_bins.at(I)});
            bin_edges.push_back(spans.at(I * 2));
            bin_edges.push_back(spans.at(I * 2 + 1));
        } else {
            const auto &bin_edges_loc = ax_bin_edges.at(I);
            axes_data.push_back({bin_edges.size(),
                                 bin_edges.size() + bin_edges_loc.size() - 1});
            bin_edges.insert(bin_edges.end(), bin_edges_loc.begin(),
                             bin_edges_loc.end());
        }
    }

    /// Call axis init for every dimension
    template <typename binnings, std::size_t... I>
    auto unroll_axis_init(
        const std::vector<scalar_type> &spans,
        const std::vector<std::size_t> &n_bins,
        const std::vector<std::vector<scalar_type>> &ax_bin_edges,
        vector_type<dindex_range> &axes_data,
        vector_type<scalar_type> &bin_edges,
        std::index_sequence<I...> /*ids*/) const {
        (axis_init<I, binnings>(spans, n_bins, ax_bin_edges, axes_data,
                                bin_edges),
         ...);
    }

    vecmem::memory_resource *m_resource{};
};

template <typename value_t, template <std::size_t> class serializer_t,
          typename populator_impl_t, typename algebra_t>
template <typename grid_t>
auto grid_factory<value_t, serializer_t, populator_impl_t,
                  algebra_t>::to_string(const grid_t &gr) const noexcept
    -> void {

    using entry_t = typename grid_t::value_type;

    // Loop over the first dimension
    const auto &ax0 = gr.template get_axis<0>();
    std::cout << "{";
    for (std::size_t i{0}; i < ax0.nbins(); ++i) {

        // Loop over the second dimension
        if constexpr (grid_t::Dim > 1) {
            const auto &ax1 = gr.template get_axis<1>();
            std::cout << "{";
            for (std::size_t j{0}; j < ax1.nbins(); ++j) {

                // Loop over the third dimension
                if constexpr (grid_t::Dim > 2) {
                    const auto &ax2 = gr.template get_axis<2>();
                    std::cout << "{";
                    for (std::size_t k{0}; k < ax2.nbins(); ++k) {

                        // Print the bin content - three dimensions
                        std::cout << "( ";
                        for (const auto &entry : gr.at(i, j, k)) {
                            if (entry == detail::invalid_value<entry_t>()) {
                                std::cout << "none ";
                            } else {
                                std::cout << entry << " ";
                            }
                        }
                        std::cout << ")";
                    }
                    std::cout << "}" << std::endl << std::endl;
                } else {
                    // Print the bin content - two dimensions
                    std::cout << "( ";
                    for (const auto &entry : gr.at(i, j)) {
                        if (entry == detail::invalid_value<entry_t>()) {
                            std::cout << "none ";
                        } else {
                            std::cout << entry << " ";
                        }
                    }
                    std::cout << ")";
                }
            }
            std::cout << "}" << std::endl;
        } else {
            // Print the bin content - one dimension
            std::cout << "( ";
            for (const auto &entry : gr.at(i)) {
                if (entry == detail::invalid_value<entry_t>()) {
                    std::cout << "none ";
                } else {
                    std::cout << entry << " ";
                }
            }
            std::cout << ")" << std::endl;
        }
    }
    std::cout << "}" << std::endl;
}

// Infer a grid factory type from an already completely assembled grid type
template <typename grid_t,
          template <std::size_t> class serializer_t =
              grid_t::template serializer_type,
          typename algebra_t = __plugin::transform3<detray::scalar>>
using grid_factory_type =
    grid_factory<typename grid_t::value_type, serializer_t,
                 typename grid_t::populator_impl, algebra_t>;

}  // namespace detray