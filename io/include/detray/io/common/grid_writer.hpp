/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/io/common/detail/definitions.hpp"
#include "detray/io/common/detail/utils.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/surface_finders/accelerator_grid.hpp"

// System include(s)
#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace detray {

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

        header_data.common = base_type::serialize(det_name, tag);

        header_data.sub_header.emplace();
        auto& grid_sub_header = header_data.sub_header.value();
        grid_sub_header.n_grids = get_n_grids(det.surface_store());

        return header_data;
    }

    /// Serialize the grid collections of a detector @param det into their io
    /// payload
    static detector_grids_payload serialize(
        const detector_t& det, const typename detector_t::name_map&) {

        detector_grids_payload grids_data;

        for (const auto& vol_desc : det.volumes()) {
            // Links to all acceleration data structures in the volume
            const auto& multi_link = vol_desc.full_link();

            for (dindex i = 0u; i < multi_link.size(); ++i) {
                const auto& acc_link = multi_link[i];
                // Don't look at empty links
                if (acc_link.is_invalid()) {
                    continue;
                }

                // If the accelerator is a grid, insert the payload
                det.surface_store().template visit<get_grid_payload>(
                    acc_link, vol_desc.index(), grids_data);
            }
        }

        return grids_data;
    }

    /// Serialize a grid @param gr of type @param type and index @param idx
    /// into its io payload
    template <class grid_t>
    static grid_payload serialize(std::size_t volume_index,
                                  io::detail::acc_type type,
                                  const std::size_t idx, const grid_t& gr) {
        grid_payload grid_data;

        grid_data.volume_link = base_type::serialize(volume_index);
        grid_data.acc_link = base_type::serialize(type, idx);

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
    /// Retrieve a @c grid_payload from grid collection elements
    struct get_grid_payload {

        template <typename grid_group_t, typename index_t>
        inline void operator()(
            [[maybe_unused]] const grid_group_t& coll,
            [[maybe_unused]] const index_t& index,
            [[maybe_unused]] std::size_t volume_index,
            [[maybe_unused]] detector_grids_payload& grids_data) const {
            using accel_t = typename grid_group_t::value_type;

            if constexpr (detail::is_grid_v<accel_t>) {

                grids_data.grids.push_back(
                    serialize(volume_index, io::detail::get_grid_id<accel_t>(),
                              index, coll[index]));
            }
        }
    };

    /// Retrieve number of overall grids in detector
    template <std::size_t I = 0u>
    static std::size_t get_n_grids(
        const typename detector_t::surface_container& store,
        std::size_t n = 0u) {

        using store_t = typename detector_t::surface_container;
        constexpr auto coll_id{store_t::value_types::to_id(I)};
        using accel_t = typename store_t::template get_type<coll_id>;

        if constexpr (detail::is_grid_v<accel_t>) {
            n += store.template size<coll_id>();
        }

        if constexpr (I < store_t::n_collections() - 1u) {
            return get_n_grids<I + 1>(store, n);
        }
        return n;
    }
};

}  // namespace detray
