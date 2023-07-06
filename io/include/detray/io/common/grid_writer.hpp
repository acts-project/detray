/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/coordinates.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/io/common/detail/utils.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/surface_finders/accelerator_grid.hpp"

// System include(s)
#include <array>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace detray {

namespace detail {

template <class grid_t>
struct is_grid : public std::false_type {};

template <typename multi_axis_t, typename value_t,
          template <std::size_t> class serializer_t, typename populator_impl_t>
struct is_grid<grid<multi_axis_t, value_t, serializer_t, populator_impl_t>>
    : public std::true_type {};

template <typename T>
inline constexpr bool is_grid_v = is_grid<T>::value;

}  // namespace detail

/// @brief Abstract base class for accelerator grid writers
template <class detector_t>
class grid_writer : public writer_interface<detector_t> {

    using base_type = writer_interface<detector_t>;

    protected:
    /// Tag the writer as "grids"
    inline static const std::string tag = "grids";

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    protected:
    /// Serialize the header information into its payload
    static grid_header_payload write_header(const detector_t& det,
                                            const std::string_view det_name) {
        grid_header_payload header_data;

        header_data.version = detail::get_detray_version();
        header_data.detector = det_name;
        header_data.tag = tag;
        header_data.date = detail::get_current_date();

        const auto& accel_store = det.surface_store();
        header_data.n_grids =
            accel_store.total_size() -
            accel_store
                .template size<detector_t::sf_finders::id::e_brute_force>();

        return header_data;
    }

    /// Serialize the grid collections of a detector @param det into their io
    /// payload
    static detector_grids_payload serialize(const detector_t& det) {

        detector_grids_payload grids_data;

        // Access the acceleration data structures recursively
        get_grid_payload(det.surface_store(), grids_data);

        return grids_data;
    }

    /// Serialize a grid @param gr of type @param type and index @param idx
    /// into its io payload
    template <class grid_t>
    static grid_payload serialize(io::detail::acc_type type,
                                  const std::size_t idx, const grid_t& gr) {
        grid_payload grid_data;

        grid_data.type = type;
        grid_data.index = idx;

        // Serialize the multi-axis into single axis payloads
        const std::array<axis_payload, grid_t::Dim> axes_data =
            serialize(gr.axes());

        grid_data.axes.resize(axes_data.size());
        std::copy(std::cbegin(axes_data), std::cend(axes_data),
                  std::begin(grid_data.axes));

        // Write the surface indices
        for (unsigned int gid = 0u; gid < gr.nbins(); ++gid) {
            // Get the local bin indices and serialize the bin into its payload
            grid_bin_payload binp = serialize(gr.serialize(gid), gr.bin(gid));
            grid_data.bins.push_back(std::move(binp));
        }

        return grid_data;
    }

    /// Serialize a multi-axis @param axes into its io payload
    template <bool ownership, typename local_frame_t, typename... axis_ts>
    static auto serialize(
        const n_axis::multi_axis<ownership, local_frame_t, axis_ts...>& axes) {

        // Serialize every single axis and construct array from their payloads
        std::array<axis_payload, sizeof...(axis_ts)> axes_data{
            serialize(axes.template get_axis<axis_ts>())...};

        return axes_data;
    }

    /// Serialize a single axis @param axis into its io payload
    template <typename bounds_t, typename binning_t>
    static axis_payload serialize(
        const n_axis::single_axis<bounds_t, binning_t>& axis) {
        axis_payload axis_data;

        axis_data.binning = axis.binning();
        axis_data.bounds = axis.bounds();
        axis_data.label = axis.label();
        axis_data.bins = axis.nbins();

        if (axis.binning() == n_axis::binning::e_regular) {
            axis_data.edges = {axis.min(), axis.max()};
        } else {
            const auto& bin_edges = axis.bin_edges();
            axis_data.edges.resize(bin_edges.size());
            std::copy(std::cbegin(bin_edges), std::cend(bin_edges),
                      std::begin(axis_data.edges));
        }

        return axis_data;
    }

    /// Serialize a multi-bin @param mbin into its io payload
    template <std::size_t DIM, typename content_range_t>
    static grid_bin_payload serialize(const n_axis::multi_bin<DIM> mbin,
                                      const content_range_t& content) {
        grid_bin_payload bin_data;

        // Local bin indices are written in the order the grid axis are stored
        for (unsigned int i = 0u; i < DIM; ++i) {
            bin_data.loc_index.push_back(mbin[i]);
        }

        // Put all entries of the bin into the payload
        bin_data.content.reserve(content.size());
        for (const auto& entry : content) {
            bin_data.content.push_back(entry.index());
        }

        return bin_data;
    }

    private:
    /// Retrieve @c grid_payload s from grid collection elements
    template <std::size_t I = 0>
    static void get_grid_payload(
        const typename detector_t::surface_container& store,
        detector_grids_payload& grids_data) {

        using store_t = typename detector_t::surface_container;
        constexpr auto coll_id{store_t::value_types::to_id(I)};
        using accelerator_t = typename store_t::template get_type<coll_id>;

        if constexpr (detail::is_grid_v<accelerator_t>) {

            const auto& coll = store.template get<coll_id>();

            for (unsigned int i = 0u; i < coll.size(); ++i) {
                grids_data.grids.push_back(serialize(
                    get_id(accelerator_t::local_frame()), i, coll[i]));
            }
        }

        if constexpr (I < store_t::n_collections() - 1u) {
            get_grid_payload<I + 1>(store, grids_data);
        }
    }

    // Infer the grid id from its coordinate system
    template <template <typename> class frame_t, typename algebra_t>
    static constexpr io::detail::acc_type get_id(const frame_t<algebra_t>&) {
        // Local coordinate frame
        using frame = frame_t<algebra_t>;

        if constexpr (std::is_same_v<frame, cartesian2<algebra_t>>) {
            return io::detail::acc_type::cartesian2_grid;
        } else if constexpr (std::is_same_v<frame, cartesian3<algebra_t>>) {
            return io::detail::acc_type::cuboid3_grid;
        } else if constexpr (std::is_same_v<frame, polar2<algebra_t>>) {
            return io::detail::acc_type::polar2_grid;
        } else if constexpr (std::is_same_v<frame, cylindrical2<algebra_t>>) {
            return io::detail::acc_type::cylinder2_grid;
        } else if constexpr (std::is_same_v<frame, cylindrical3<algebra_t>>) {
            return io::detail::acc_type::cylinder3_grid;
        } else {
            return io::detail::acc_type::unknown;
        }
    }
};

}  // namespace detray
