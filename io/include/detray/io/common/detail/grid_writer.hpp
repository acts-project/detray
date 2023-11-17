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
#include "detray/utils/type_list.hpp"

// System include(s)
#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace detray::detail {

/// @brief Abstract base class for accelerator grid writers
template <typename detector_t, typename value_t>
class grid_writer : public writer_interface<detector_t> {

    using base_type = writer_interface<detector_t>;

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    /// Serialize the header information into its payload
    template <typename grid_store_t>
    static grid_header_payload write_header(const std::string_view writer_tag,
                                            const grid_store_t& store,
                                            const std::string_view det_name) {

        grid_header_payload header_data;

        header_data.common = base_type::serialize(det_name, writer_tag);

        header_data.sub_header.emplace();
        auto& grid_sub_header = header_data.sub_header.value();
        grid_sub_header.n_grids = get_n_grids(store);

        return header_data;
    }

    protected:
    /// Serialize a grid @param gr of type @param type and index @param idx
    /// into its io payload, using @param serializer for the bin content
    template <typename content_t, class grid_t>
    static grid_payload<content_t> serialize(
        std::size_t volume_index, io::detail::acc_type type,
        const std::size_t idx, const grid_t& gr,
        std::function<content_t(const value_t&)> serializer) {

        grid_payload<content_t> grid_data;

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
            grid_bin_payload binp =
                serialize(gr.serialize(gid), gr.bin(gid), serializer);
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
    template <typename content_t, std::size_t DIM, typename content_range_t>
    static grid_bin_payload<content_t> serialize(
        const n_axis::multi_bin<DIM> mbin, const content_range_t& content,
        std::function<content_t(const value_t&)> serializer) {

        grid_bin_payload<content_t> bin_data;

        // Local bin indices are written in the order the grid axis are stored
        for (unsigned int i = 0u; i < DIM; ++i) {
            bin_data.loc_index.push_back(mbin[i]);
        }

        // Put all entries of the bin into the payload
        bin_data.content.reserve(content.size());
        for (const auto& entry : content) {
            bin_data.content.push_back(serializer(entry));
        }

        return bin_data;
    }

    /// Serialize a grid from a collection into its payload
    ///
    /// @param store the data store of grids (tuple of grid collections)
    /// @param grid_link type and index of the grid
    /// @param owner_idx inder of the owner of the grid (e.g. volume index)
    /// @param grid_data the grid payload to be filled
    /// @param serializer callable that can serialize a grid bin entry into its
    /// respective IO payload (of type @tparam content_t)
    template <typename store_t, typename content_t, typename serializer_t>
    static void serialize(const store_t& store,
                          typename store_t::single_link grid_link,
                          dindex owner_idx,
                          detector_grids_payload<content_t>& grids_data,
                          serializer_t serializer) {

        // If the accelerator is a grid, insert the payload
        store.template visit<get_grid_payload>(grid_link, owner_idx, grids_data,
                                               serializer);
    }

    private:
    /// Retrieve a @c grid_payload from grid collection elements
    struct get_grid_payload {

        template <typename grid_group_t, typename index_t, typename content_t,
                  typename serializer_t>
        inline void operator()(
            [[maybe_unused]] const grid_group_t& coll,
            [[maybe_unused]] const index_t& index,
            [[maybe_unused]] std::size_t link,
            [[maybe_unused]] detector_grids_payload<content_t>& grids_data,
            [[maybe_unused]] serializer_t& serializer) const {

            using accel_t = typename grid_group_t::value_type;

            if constexpr (detail::is_grid_v<accel_t>) {

                grids_data.grids.push_back(serialize<content_t>(
                    link, io::detail::get_grid_id<accel_t>(), index,
                    coll[index], serializer));
            }
        }
    };

    /// Retrieve number of overall grids in detector
    template <std::size_t I = 0u, typename store_t>
    static std::size_t get_n_grids(const store_t& store, std::size_t n = 0u) {

        constexpr auto coll_id{store_t::value_types::to_id(I)};
        using value_type = typename store_t::template get_type<coll_id>;

        if constexpr (detail::is_grid_v<value_type>) {
            n += store.template size<coll_id>();
        }

        if constexpr (I < store_t::n_collections() - 1u) {
            return get_n_grids<I + 1>(store, n);
        }
        return n;
    }
};

}  // namespace detray::detail
