/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/io/common/detail/type_traits.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/tools/grid_builder.hpp"
#include "detray/tools/grid_factory.hpp"
#include "detray/utils/ranges.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <queue>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace detray {

/// @brief Abstract base class for surface grid readers
template <class detector_t,
          typename value_t = typename detector_t::surface_type,
          typename CAP = std::integral_constant<std::size_t, 9>,
          typename DIM = std::integral_constant<std::size_t, 2>>
class grid_reader : public reader_interface<detector_t> {

    using base_type = reader_interface<detector_t>;
    /// IO accelerator ids do not need to coincide with the detector ids,
    /// because they are shared with ACTS
    using acc_type = io::detail::acc_type;
    using algebra_t = typename detector_t::transform3;
    using scalar_t = typename algebra_t::scalar_type;

    // Bin content type to be written into the grid bin payloads
    // For detector surface descriptors, write only the surface index
    using content_t = std::conditional_t<
        std::is_same_v<value_t, typename detector_t::surface_type>, std::size_t,
        value_t>;

    static constexpr std::size_t dim{DIM()};
    static constexpr std::size_t bin_capacity{CAP()};

    protected:
    /// Tag the reader as "surface_grids"
    inline static const std::string tag = "surface_grids";

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    protected:
    /// Deserialize the detector grids @param grids_data from their IO
    /// payload
    static void deserialize(
        detector_builder<typename detector_t::metadata, volume_builder>
            &det_builder,
        typename detector_t::name_map &,
        const detector_grids_payload<content_t> &grids_data) {

        // Deserialize the grids volume by volume
        for (const auto &grid_data : grids_data.grids) {

            std::queue<axis::bounds> bounds;
            std::queue<axis::binning> binnings;
            for (const auto &axis_data : grid_data.axes) {
                bounds.push(axis_data.bounds);
                binnings.push(axis_data.binning);
            }
            deserialize(bounds, binnings, grid_data, det_builder);
        }
    }

    private:
    /// @brief recursively build the grid: axis bounds (open, closed, circular)
    ///
    /// @tparam bounds_ts type list that contains the bounds types that were
    ///         identified from the IO ids so far (start with empty list)
    /// @tparam binning_ts type list that contains the binning types that were
    ///         identified from the IO ids so far (start with empty list)
    ///
    /// @param bound_ids runtime queue of bounds type ids (read from file)
    /// @param binning_ids runtime queue of binning type ids (read from file)
    template <typename bounds_ts = types::list<>,
              typename binning_ts = types::list<>, typename... Ts>
    static void deserialize(std::queue<axis::bounds> &bound_ids,
                            std::queue<axis::binning> &binning_ids,
                            Ts &&... data) {
        using namespace axis;

        constexpr std::size_t n_bounds_types{types::size<bounds_ts>};

        // Base case: If the bounds types are filled, continue with the binnings
        if constexpr (n_bounds_types == dim) {
            return deserialize<bounds_ts, binning_ts>(
                binning_ids, std::forward<Ts>(data)...);
        } else if (!bound_ids.empty()) {
            // The axis label, e.g. x, y or z by number
            constexpr auto lb{static_cast<label>(n_bounds_types)};

            const auto first_id{bound_ids.front()};
            bound_ids.pop();

            // Based on the type id, add the next bounds type to the type list
            // and continue
            switch (first_id) {
                case bounds::e_closed: {
                    using new_bounds_ts =
                        types::push_back<bounds_ts, closed<lb>>;
                    return deserialize<new_bounds_ts, binning_ts>(
                        bound_ids, binning_ids, std::forward<Ts>(data)...);
                }
                case bounds::e_open: {
                    using new_bounds_ts = types::push_back<bounds_ts, open<lb>>;
                    return deserialize<new_bounds_ts, binning_ts>(
                        bound_ids, binning_ids, std::forward<Ts>(data)...);
                }
                case bounds::e_circular: {
                    using new_bounds_ts =
                        types::push_back<bounds_ts, circular<lb>>;
                    return deserialize<new_bounds_ts, binning_ts>(
                        bound_ids, binning_ids, std::forward<Ts>(data)...);
                }
                // Test some edge cases
                default: {
                    throw std::invalid_argument(
                        "Given type id could not be matched to a grid boundary "
                        "type: " +
                        std::to_string(static_cast<std::int64_t>(first_id)));
                    break;
                }
            };
        }
    }

    /// @brief recursively build the grid: axis binning (regular, irregular)
    ///
    /// @tparam bounds_ts type list that contains the bounds types
    /// @tparam binning_ts type list that contains the binning types that were
    ///         identified from the IO ids so far (start with empty list)
    ///
    /// @param binning_ids runtime queue of binning type ids (read from file)
    template <typename bounds_ts, typename binning_ts, typename... Ts,
              std::enable_if_t<types::size<bounds_ts> == dim, bool> = true>
    static void deserialize(std::queue<axis::binning> &binning_ids,
                            Ts &&... data) {

        using namespace axis;

        using regular_binning_t = regular<host_container_types, scalar_t>;
        using irregular_binning_t = irregular<host_container_types, scalar_t>;

        // Base case: If the binning types are filled, continue with the frame
        if constexpr (types::size<binning_ts> == dim) {
            return deserialize<bounds_ts, binning_ts>(
                std::forward<Ts>(data)...);
        } else if (!binning_ids.empty()) {

            const auto first_id{binning_ids.front()};
            binning_ids.pop();

            switch (first_id) {
                case binning::e_regular: {
                    using new_binning_ts =
                        types::push_back<binning_ts, regular_binning_t>;
                    return deserialize<bounds_ts, new_binning_ts>(
                        binning_ids, std::forward<Ts>(data)...);
                }
                case binning::e_irregular: {
                    using new_binning_ts =
                        types::push_back<binning_ts, irregular_binning_t>;
                    return deserialize<bounds_ts, new_binning_ts>(
                        binning_ids, std::forward<Ts>(data)...);
                }
                // Test some edge cases
                default: {
                    throw std::invalid_argument(
                        "Given type id could not be matched to a grid binning "
                        "type: " +
                        std::to_string(static_cast<std::int64_t>(first_id)));
                    break;
                }
            };
        }
    }

    /// @brief recursively build the grid: find the grid geometry (loc. coord.)
    ///
    /// @tparam bounds_ts type list that contains the bounds types
    /// @tparam binning_ts type list that contains the binning types
    ///
    /// @param grid_data grid IO payload (read from file)
    /// @param det_builder gather the grid data and build the final volume
    template <typename bounds_ts, typename binning_ts,
              std::enable_if_t<types::size<bounds_ts> == dim and
                                   types::size<binning_ts> == dim,
                               bool> = true>
    static void deserialize(const grid_payload<content_t> &grid_data,
                            detector_builder<typename detector_t::metadata,
                                             volume_builder> &det_builder) {

        // Throw expection if the accelerator link type id is invalid
        auto print_error = [](io::detail::acc_type acc_link) -> void {
            if (acc_link == io::detail::acc_type::unknown) {
                throw std::invalid_argument(
                    "Unknown accelerator id in geometry file!");
            } else {
                throw std::invalid_argument(
                    "Given accelerator id could not be matched to a grid "
                    "type: " +
                    std::to_string(static_cast<std::int64_t>(acc_link)));
            }
        };

        // Need to pass actual instances to deduce contained types as
        // template parameter packs
        constexpr auto bounds = bounds_ts{};
        constexpr auto binnings = binning_ts{};

        // Check only 2-dimensional grid types
        if constexpr (dim >= 2) {
            switch (grid_data.acc_link.type) {
                // rectangle, trapezoid, (triangle) grids
                case io::detail::acc_type::cartesian2_grid: {
                    return deserialize<cartesian2<algebra_t>>(
                        grid_data, det_builder, bounds, binnings);
                }
                // ring/disc, annulus grids
                case io::detail::acc_type::polar2_grid: {
                    return deserialize<polar2<algebra_t>>(
                        grid_data, det_builder, bounds, binnings);
                }
                // 2D cylinder grid
                case io::detail::acc_type::cylinder2_grid: {
                    return deserialize<cylindrical2<algebra_t>>(
                        grid_data, det_builder, bounds, binnings);
                }
                default: {
                    print_error(grid_data.acc_link.type);
                    break;
                }
            };
        } else if constexpr (dim >= 3) {
            switch (grid_data.acc_link.type) {
                // cuboid grid
                case io::detail::acc_type::cuboid3_grid: {
                    return deserialize<cartesian3<algebra_t>>(
                        grid_data, det_builder, bounds, binnings);
                }
                // 3D cylinder grid
                case io::detail::acc_type::cylinder3_grid: {
                    return deserialize<cylindrical3<algebra_t>>(
                        grid_data, det_builder, bounds, binnings);
                }
                default: {
                    print_error(grid_data.acc_link.type);
                    break;
                }
            };
        } else {
            throw std::invalid_argument("No 1D grid type defined in detray");
        }
    }

    /// @brief End of recursion: build the grid from the @param grid_data
    template <typename local_frame_t, typename... bounds_ts,
              typename... binning_ts,
              std::enable_if_t<sizeof...(bounds_ts) == dim and
                                   sizeof...(binning_ts) == dim,
                               bool> = true>
    static void deserialize(const grid_payload<content_t> &grid_data,
                            detector_builder<typename detector_t::metadata,
                                             volume_builder> &det_builder,
                            types::list<bounds_ts...>,
                            types::list<binning_ts...>) {
        // Assemble the grid type
        using axes_t =
            axis::multi_axis<false, local_frame_t,
                             axis::single_axis<bounds_ts, binning_ts>...>;
        // Ok for now: Only have this serializer
        using grid_t = grid<axes_t, bins::static_array<value_t, bin_capacity>,
                            simple_serializer>;

        static_assert(grid_t::dim == dim,
                      "Grid dimension does not meet dimension of grid reader");

        const auto volume_idx{base_type::deserialize(grid_data.volume_link)};

        // Error output
        std::stringstream err_stream;
        err_stream << "Volume " << volume_idx << ": ";

        // The compiler will instantiate this function for all possible types of
        // grids: Only proceed, if the grid type is known by the detector
        if constexpr (detector_t::accel::template is_defined<grid_t>()) {

            // Decorate the current volume builder with the grid
            using grid_builder_t =
                grid_builder<detector_t, grid_t, detail::fill_by_pos>;
            auto v_builder =
                det_builder.template decorate<grid_builder_t>(volume_idx);
            auto vgr_builder = dynamic_cast<grid_builder_t *>(v_builder);

            // Initialize the grid axes
            std::vector<std::size_t> n_bins_per_axis{};
            std::vector<scalar_t> spans{};
            std::vector<std::vector<scalar_t>> ax_bin_edges{};

            for (const auto &axis_data : grid_data.axes) {
                n_bins_per_axis.push_back(axis_data.bins);
                std::vector<scalar_t> edges{};
                std::copy(axis_data.edges.begin(), axis_data.edges.end(),
                          std::back_inserter(edges));
                ax_bin_edges.emplace_back(std::move(edges));
                spans.push_back(static_cast<scalar_t>(axis_data.edges.front()));
                spans.push_back(static_cast<scalar_t>(axis_data.edges.back()));
            }

            vgr_builder->init_grid(spans, n_bins_per_axis, {}, ax_bin_edges);
            auto &grid = vgr_builder->get();
            const std::size_t n_bins{grid.nbins()};

            value_t empty_sf{};
            axis::multi_bin<dim> mbin;
            for (const auto &bin_data : grid_data.bins) {

                assert(dim == bin_data.loc_index.size() &&
                       "Numer of local bin indices in input file does not "
                       "match grid dimension");

                // The local bin indices for the bin to be filled
                for (const auto &[i, bin_idx] :
                     detray::views::enumerate(bin_data.loc_index)) {
                    mbin[i] = bin_idx;
                }

                const auto gbin = grid.serializer()(grid.axes(), mbin);
                if (gbin >= n_bins) {
                    err_stream << "Bin index " << mbin << " out of bounds";
                    throw std::invalid_argument(err_stream.str());
                }

                // For now assume surfaces ids as the only grid input
                for (const auto c : bin_data.content) {
                    empty_sf.set_volume(volume_idx);
                    empty_sf.set_index(static_cast<dindex>(c));
                    vgr_builder->get().template populate<attach<>>(mbin,
                                                                   empty_sf);
                }
            }
        } else {
            types::print<types::list<grid_t>>();
            err_stream
                << "Grid type in file does not match any grid type in detector";
            throw std::invalid_argument(err_stream.str());
        }
    }
};

}  // namespace detray
