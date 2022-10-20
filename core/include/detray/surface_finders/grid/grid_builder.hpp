/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/containers.hpp"
#include "detray/masks/masks.hpp"
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/grid.hpp"
#include "detray/surface_finders/grid/grid_collection.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <iostream>
#include <string>

namespace detray {

namespace detail {

/// @brief Helper type to assemble an multi-axis from shapes tuple and binnings
template <bool is_owning, typename containers, typename local_frame, typename,
          typename>
struct multi_axis_assembler;

/// @brief Specialized struct to extract axis shapes from a tuple
template <bool is_owning, typename containers, typename local_frame,
          typename... axis_shapes, typename... binning_ts>
struct multi_axis_assembler<is_owning, containers, local_frame,
                            dtuple<axis_shapes...>, dtuple<binning_ts...>> {

    static_assert(sizeof...(axis_shapes) == sizeof...(binning_ts),
                  "Number of axis shapes for this mask and given binning types "
                  "don't match!");

    using type =
        n_axis::multi_axis<is_owning, local_frame,
                           n_axis::single_axis<axis_shapes, binning_ts>...>;
};

}  // namespace detail

/// Typedef for easier construction from mask shapes
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
          typename algebra_t = __plugin::transform3<detray::scalar>,
          typename container_t = host_container_types>
class grid_builder {

    public:
    // All grids are owning since they are used to fill the data
    static constexpr bool is_owning = true;

    using value_type = value_t;
    template <typename grid_shape_t>
    using grid_type =
        grid<grid_shape_t, value_type, serializer_t, populator_impl_t>;

    using scalar_type = typename algebra_t::scalar_type;
    template <typename T>
    using vector_type = typename container_t::template vector_type<T>;

    /// Takes the resource of the detector to allocate memory correctly
    grid_builder(vecmem::memory_resource &resource) : m_resource(&resource) {}

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

    /// @brief Build empty grids from boundary values.
    struct grid_factory {

        //
        // annulus 2D
        //
        template <
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_r = n_axis::regular,
            template <typename, typename> class binning_phi = n_axis::regular>
        auto new_grid(
            const mask<annulus2D<>> &grid_bounds, const std::size_t n_bins_r,
            const std::size_t n_bins_phi,
            const std::vector<scalar_type> &bin_edges_r = {},
            const std::vector<scalar_type> &bin_edges_phi = {}) const {
            // Axes boundaries
            using boundary = annulus2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<annulus2D<>, e_shape, binning_r, binning_phi>(
                {b_values[boundary::e_min_r], b_values[boundary::e_max_r],
                 b_values[boundary::e_average_phi] -
                     b_values[boundary::e_min_phi_rel],
                 b_values[boundary::e_average_phi] +
                     b_values[boundary::e_max_phi_rel]},
                {n_bins_r, n_bins_phi}, {bin_edges_r, bin_edges_phi});
        }

        //
        // cuboid 3D
        //
        template <
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_x = n_axis::regular,
            template <typename, typename> class binning_y = n_axis::regular,
            template <typename, typename> class binning_z = n_axis::regular>
        auto new_grid(const mask<cuboid3D> &grid_bounds,
                      const std::size_t n_bins_x, const std::size_t n_bins_y,
                      const std::size_t n_bins_z,
                      const std::vector<scalar_type> &bin_edges_x = {},
                      const std::vector<scalar_type> &bin_edges_y = {},
                      const std::vector<scalar_type> &bin_edges_z = {}) const {
            // Axes boundaries
            using boundary = cuboid3D::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<cuboid3D, e_shape, binning_x, binning_y, binning_z>(
                {-b_values[boundary::e_half_x], b_values[boundary::e_half_x],
                 -b_values[boundary::e_half_y], b_values[boundary::e_half_y],
                 -b_values[boundary::e_half_z], b_values[boundary::e_half_z]},
                {n_bins_x, n_bins_y, n_bins_z},
                {bin_edges_x, bin_edges_y, bin_edges_z});
        }

        //
        // cylinder 2D
        //
        template <
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_phi = n_axis::regular,
            template <typename, typename> class binning_z = n_axis::regular>
        auto new_grid(const mask<cylinder2D<>> &grid_bounds,
                      const std::size_t n_bins_phi, const std::size_t n_bins_z,
                      const std::vector<scalar_type> &bin_edges_phi = {},
                      const std::vector<scalar_type> &bin_edges_z = {}) const {
            // Axes boundaries
            using boundary = cylinder2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<cylinder2D<>, e_shape, binning_phi, binning_z>(
                {0.f, 2.f * static_cast<scalar_type>(M_PI),
                 b_values[boundary::e_n_half_z],
                 b_values[boundary::e_p_half_z]},
                {n_bins_phi, n_bins_z}, {bin_edges_phi, bin_edges_z});
        }

        //
        // cylinder 3D
        //
        template <
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_r = n_axis::regular,
            template <typename, typename> class binning_phi = n_axis::regular,
            template <typename, typename> class binning_z = n_axis::regular>
        auto new_grid(const mask<cylinder3D> &grid_bounds,
                      const std::size_t n_bins_r, const std::size_t n_bins_phi,
                      const std::size_t n_bins_z,
                      const std::vector<scalar_type> &bin_edges_r = {},
                      const std::vector<scalar_type> &bin_edges_phi = {},
                      const std::vector<scalar_type> &bin_edges_z = {}) const {
            // Axes boundaries
            using boundary = cylinder3D::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<cylinder3D, e_shape, binning_r, binning_phi,
                            binning_z>(
                {0.f, b_values[boundary::e_r], 0.f,
                 2.f * static_cast<scalar_type>(M_PI),
                 -b_values[boundary::e_n_half_z],
                 b_values[boundary::e_p_half_z]},
                {n_bins_r, n_bins_phi, n_bins_z},
                {bin_edges_r, bin_edges_phi, bin_edges_z});
        }

        //
        // polar 2D
        //
        template <
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_phi = n_axis::regular,
            template <typename, typename> class binning_z = n_axis::regular>
        auto new_grid(
            const mask<ring2D<>> &grid_bounds, const std::size_t n_bins_r,
            const std::size_t n_bins_phi,
            const std::vector<scalar_type> &bin_edges_r = {},
            const std::vector<scalar_type> &bin_edges_phi = {}) const {
            // Axes boundaries
            using boundary = ring2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<ring2D<>, e_shape, binning_phi, binning_z>(
                {b_values[boundary::e_inner_r], b_values[boundary::e_outer_r],
                 0.f, 2.f * static_cast<scalar_type>(M_PI)},
                {n_bins_r, n_bins_phi}, {bin_edges_r, bin_edges_phi});
        }

        //
        // rectangle 2D
        //
        template <
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_phi = n_axis::regular,
            template <typename, typename> class binning_z = n_axis::regular>
        auto new_grid(const mask<rectangle2D<>> &grid_bounds,
                      const std::size_t n_bins_x, const std::size_t n_bins_y,
                      const std::vector<scalar_type> &bin_edges_x = {},
                      const std::vector<scalar_type> &bin_edges_y = {}) const {
            // Axes boundaries
            using boundary = rectangle2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<rectangle2D<>, e_shape, binning_phi, binning_z>(
                {-b_values[boundary::e_half_x], b_values[boundary::e_half_x],
                 -b_values[boundary::e_half_y], b_values[boundary::e_half_y]},
                {n_bins_x, n_bins_y}, {bin_edges_x, bin_edges_y});
        }

        //
        // trapezoid 2D
        //
        template <
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_phi = n_axis::regular,
            template <typename, typename> class binning_z = n_axis::regular>
        auto new_grid(const mask<trapezoid2D<>> &grid_bounds,
                      const std::size_t n_bins_x, const std::size_t n_bins_y,
                      const std::vector<scalar_type> &bin_edges_x = {},
                      const std::vector<scalar_type> &bin_edges_y = {}) const {
            // Axes boundaries
            using boundary = trapezoid2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<trapezoid2D<>, e_shape, binning_phi, binning_z>(
                {-b_values[boundary::e_half_length_1],
                 b_values[boundary::e_half_length_1],
                 -b_values[boundary::e_half_length_2],
                 b_values[boundary::e_half_length_2]},
                {n_bins_x, n_bins_y}, {bin_edges_x, bin_edges_y});
        }

        /// @brief Create and empty grid with fully initialized axes.
        ///
        /// @tparam grid_shape_t the shape of the resulting grid
        ///         (e.g. cylinder2D).
        /// @tparam e_shape the shape of the regular axes
        ///         (open vs. closed binning).
        /// @tparam binnings the binning types of the axes
        ///         (regular vs. irregular)
        ///
        /// @param spans the span of the axis values for regular axes, otherwise
        ///              ignored.
        /// @param n_bins the number of bins for regular axes, otherwise ignored
        /// @param ax_bin_edges the explicit bin edges for irregular axes
        ///                     (lower bin edges + the the upper edge of the
        ///                     last bin), otherwise ignored.
        template <typename grid_shape_t,
                  n_axis::shape e_shape = n_axis::shape::e_open,
                  template <typename, typename> class... binning_ts>
        auto new_grid(const std::vector<scalar_type> spans,
                      const std::vector<std::size_t> n_bins,
                      const std::vector<std::vector<scalar_type>>
                          &ax_bin_edges = {{}}) const {

            // Build the coordinate axes and the grid
            using axes_t = coordinate_axes<
                typename grid_shape_t::template axes<e_shape, binning_ts...>,
                is_owning, container_t, algebra_t>;
            using bin_t = typename grid_type<axes_t>::bin_type;

            /// make binning types accessible
            using binnings =
                std::tuple<binning_ts<host_container_types, scalar_type>...>;

            // Prepare data
            vector_type<dindex_range> axes_data(m_resource);
            vector_type<scalar_type> bin_edges(m_resource);

            // Call init for every axis
            unroll_axis_init<binnings>(spans, n_bins, ax_bin_edges, axes_data,
                                       bin_edges,
                                       std::make_index_sequence<axes_t::Dim>{});

            // Assemble the grid and return it
            axes_t axes(std::move(axes_data), std::move(bin_edges));

            vector_type<bin_t> bin_data(m_resource);
            bin_data.resize(axes.nbins()[0] * axes.nbins()[1],
                            populator_impl_t::template init<
                                typename grid_type<axes_t>::value_type>());

            return grid_type<axes_t>(std::move(bin_data), std::move(axes));
        }

        /// Initialize a single axis (either regular or irregular)
        /// @note change to template lambda as soon as it becomes available.
        template <std::size_t I, typename binnings>
        auto axis_init(
            [[maybe_unused]] const std::vector<scalar_type> &spans,
            [[maybe_unused]] const std::vector<std::size_t> &n_bins,
            [[maybe_unused]] const std::vector<std::vector<scalar_type>>
                &ax_bin_edges,
            vector_type<dindex_range> &axes_data,
            vector_type<scalar_type> &bin_edges) const {
            if constexpr (std::is_same_v<std::tuple_element_t<I, binnings>,
                                         n_axis::regular<>>) {
                axes_data.push_back({bin_edges.size(), n_bins.at(I)});
                bin_edges.push_back(spans.at(I * 2));
                bin_edges.push_back(spans.at(I * 2 + 1));
            } else {
                const auto &bin_edges_loc = ax_bin_edges.at(I);
                axes_data.push_back(
                    {bin_edges.size(),
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

        /// The memory resource with which to allocate the grid memory
        vecmem::memory_resource *m_resource{nullptr};
    };

    /// Get new grid builder
    constexpr auto new_factory() const -> grid_factory { return {m_resource}; }

    vecmem::memory_resource *m_resource{};
};

template <typename value_t, template <std::size_t> class serializer_t,
          typename populator_impl_t, typename algebra_t, typename container_t>
template <typename grid_t>
auto grid_builder<value_t, serializer_t, populator_impl_t, algebra_t,
                  container_t>::to_string(const grid_t &gr) const noexcept
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
                                std::cout << "inv ";
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
                            std::cout << "inv ";
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
                    std::cout << "inv ";
                } else {
                    std::cout << entry << " ";
                }
            }
            std::cout << ")" << std::endl;
        }
    }
    std::cout << "}" << std::endl;
}

}  // namespace detray