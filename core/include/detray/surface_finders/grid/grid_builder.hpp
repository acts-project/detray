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

/// A two-dimensional grid for object storage
///
/// @tparam populator_t  is a prescription what to do when a bin gets
/// pupulated, it broadcasts also the value type
/// @tparam serialzier_t  type of the serializer to the storage represenations
template <typename value_t, template <std::size_t> class serializer_t,
          typename populator_impl_t,
          typename algebra_t = __plugin::transform3<detray::scalar>,
          typename container_t = host_container_types>
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
    using vector_type = typename container_t::template vector_type<T>;

    /// Takes the resource of the detector to allocate memory correctly
    grid_factory(vecmem::memory_resource &resource) : m_resource(&resource) {}

    /// Get new grid collection
    /*template<typename grid_t>
    auto new_collection() -> grid_collection<grid_t> {

    }

    /// Get new grid collection
    template<typename grid_t>
    auto to_string(const grid_t & gr) -> std::string {

    }*/

    /// @brief Build empty grids from boundary values.
    struct grid_builder {

        //
        // annulus 2D
        //
        template <
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_r = n_axis::regular,
            template <typename, typename> class binning_phi = n_axis::regular>
        auto new_grid(const mask<annulus2D<>> &grid_bounds,
                      const std::size_t n_bins_r, const std::size_t n_bins_phi,
                      const std::vector<scalar_type> &bin_edges_r = {},
                      const std::vector<scalar_type> &bin_edges_phi = {}) {
            // Axes boundaries
            using boundary = annulus2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<annulus2D<>, e_shape, binning_r, binning_phi>(
                b_values[boundary::e_min_r], b_values[boundary::e_max_r],
                b_values[boundary::e_average_phi] -
                    b_values[boundary::e_min_phi_rel],
                b_values[boundary::e_average_phi] +
                    b_values[boundary::e_max_phi_rel],
                n_bins_r, n_bins_phi, bin_edges_r, bin_edges_phi);
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
                      const std::vector<scalar_type> &bin_edges_z = {}) {
            // Axes boundaries
            using boundary = cylinder2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<cylinder2D<>, e_shape, binning_phi, binning_z>(
                0.f, 2.f * static_cast<scalar_type>(M_PI),
                b_values[boundary::e_n_half_z], b_values[boundary::e_p_half_z],
                n_bins_phi, n_bins_z, bin_edges_phi, bin_edges_z);
        }

        //
        // polar 2D
        //
        template <
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_phi = n_axis::regular,
            template <typename, typename> class binning_z = n_axis::regular>
        auto new_grid(const mask<ring2D<>> &grid_bounds,
                      const std::size_t n_bins_r, const std::size_t n_bins_phi,
                      const std::vector<scalar_type> &bin_edges_r = {},
                      const std::vector<scalar_type> &bin_edges_phi = {}) {
            // Axes boundaries
            using boundary = ring2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<ring2D<>, e_shape, binning_phi, binning_z>(
                b_values[boundary::e_inner_r], b_values[boundary::e_outer_r],
                0.f, 2.f * static_cast<scalar_type>(M_PI), n_bins_r, n_bins_phi,
                bin_edges_r, bin_edges_phi);
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
                      const std::vector<scalar_type> &bin_edges_y = {}) {
            // Axes boundaries
            using boundary = rectangle2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<rectangle2D<>, e_shape, binning_phi, binning_z>(
                -b_values[boundary::e_half_x], b_values[boundary::e_half_x],
                -b_values[boundary::e_half_y], b_values[boundary::e_half_y],
                n_bins_x, n_bins_y, bin_edges_x, bin_edges_y);
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
                      const std::vector<scalar_type> &bin_edges_y = {}) {
            // Axes boundaries
            using boundary = trapezoid2D<>::boundaries;
            auto b_values = grid_bounds.values();

            return new_grid<trapezoid2D<>, e_shape, binning_phi, binning_z>(
                -b_values[boundary::e_half_length_1],
                b_values[boundary::e_half_length_1],
                -b_values[boundary::e_half_length_2],
                b_values[boundary::e_half_length_2], n_bins_x, n_bins_y,
                bin_edges_x, bin_edges_y);
        }
        /*new_grid(const mask<cuboid3D> & grid_bounds) {

        }

        new_grid(const mask<cylinder3D> & grid_bounds) {

        }*/

        template <
            typename grid_shape_t,
            n_axis::shape e_shape = n_axis::shape::e_open,
            template <typename, typename> class binning_loc0 = n_axis::regular,
            template <typename, typename> class binning_loc1 = n_axis::regular>
        auto new_grid(const scalar min_loc0, const scalar max_loc0,
                      const scalar min_loc1, const scalar max_loc1,
                      const std::size_t n_loc0, const std::size_t n_loc1,
                      const std::vector<scalar_type> &bin_edges_loc0 = {},
                      const std::vector<scalar_type> &bin_edges_loc1 = {}) {
            // Prepare data
            vector_type<dindex_range> axes_data(m_resource);
            vector_type<scalar_type> bin_edges(m_resource);

            // Setup empty grid data

            // local 0 axis (r)
            if constexpr (std::is_same_v<
                              binning_loc0<host_container_types, scalar_type>,
                              n_axis::regular<>>) {
                axes_data.push_back({0, n_loc0});
                bin_edges.push_back(min_loc0);
                bin_edges.push_back(max_loc0);
            } else {
                axes_data.push_back({0, bin_edges_loc0.size()});
                bin_edges.insert(bin_edges.begin(), bin_edges_loc0.begin(),
                                 bin_edges_loc0.end());
            }
            // local 1 axis (phi)
            if constexpr (std::is_same_v<
                              binning_loc1<host_container_types, scalar_type>,
                              n_axis::regular<>>) {
                axes_data.push_back({bin_edges.size(), n_loc1});
                bin_edges.push_back(min_loc1);
                bin_edges.push_back(max_loc1);
            } else {
                axes_data.push_back(
                    {bin_edges.size(),
                     bin_edges.size() + bin_edges_loc1.size() - 1});
                bin_edges.insert(bin_edges.end(), bin_edges_loc1.begin(),
                                 bin_edges_loc1.end());
            }

            // Build the coordinate axes and the grid
            using axes_t =
                coordinate_axes<typename grid_shape_t::template axes<
                                    e_shape, binning_loc0, binning_loc1>,
                                is_owning, container_t, algebra_t>;

            axes_t axes(std::move(axes_data), std::move(bin_edges));

            using bin_t = typename grid_type<axes_t>::bin_type;

            vector_type<bin_t> bin_data(m_resource);
            bin_data.resize(axes.nbins()[0] * axes.nbins()[1],
                            populator_impl_t::template init<
                                typename grid_type<axes_t>::value_type>());

            return grid_type<axes_t>(std::move(bin_data), std::move(axes));
        }

        /// The memory resource with which to allocate the grid memory
        vecmem::memory_resource *m_resource{nullptr};
    };

    /// Get new grid builder
    auto new_grid_builder() -> grid_builder { return {m_resource}; }

    vecmem::memory_resource *m_resource{};
};

}  // namespace detray