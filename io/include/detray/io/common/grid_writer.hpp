/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/io/common/detail/utils.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"

// System include(s)
#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace detray {

/// @brief Abstract base class for tracking geometry writers
template <class detector_t, class grid_t>
class grid_writer : public writer_interface<detector_t> {

    using base_type = writer_interface<detector_t>;

    protected:
    /// Tag the writer as "geometry"
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

        header_data.n_grids = det.volumes().size();

        return header_data;
    }

    /// Serialize the material description of a detector @param det into its io
    /// payload
    static detector_grids_payload serialize(
        const detector_t& det) {
        using sf_finders = typename detector_t::material_container::value_types;

        detector_grids_payload grids_data;

        // The surface store does not only contain grids (can be anything)
        const auto& accel_store = det.surface_store();

        // Go through all potentially known grids types
        for (const auto [idx, gr] : detray::views::enumerate(
                 accel_store.template get<sf_finders::e_disc_grid>())) {
            grids_data.grids.push_back(serialize(gr, idx));
        }

        return grids_data;
    }

    /// Serialize a link @param idx into its io payload
    static single_link_payload serialize(const std::size_t idx) {
        single_link_payload link_data;
        link_data.link = idx;

        return link_data;
    }

    /// Serialize a detector @param det into its io payload
    static grid_payload serialize(const grid_t& gr) {
        grid_payload grid_data;

        // Serialize the multi-axis into single axis payloads
        const std::array<axis_payload, grid_t::Dim> axes_data = serialize(gr.axes());

        // Copy into payload
        grid_data.axes.resize(axes_data.size());
        std::copy(std::cbegin(axes_data), std::cend(axes_data),
                  std::begin(grid_data.axes));

        // Write the surface indices
        for (const auto& entry : gr.all()) {
            grid_data.entries.push_back(serialize(entry.index()));
        }

        return grid_data;
    }

    /// Serialize a multi-axis @param axes into its io payload
    template <bool ownership, typename local_frame_t, typename... axis_ts>
    static std::array<axis_payload, grid_t::Dim> serialize(const multi_axis<ownership, local_frame_t, axis_ts...> &axes) {

        // Serialize every single axis and construct array from their payloads
        std::array<axis_payload, grid_t::Dim> axes_data{serialize(axes.template get_axis<axis_ts>())...};

        return axes_data;
    }

    /// Serialize a single axis @param axis into its io payload
    template <typename bounds_t, typename binning_t>
    static axis_payload serialize(const single_axis<bounds_t, binning_t> &axis) {
        axis_payload axis_data;

        axis_data.binning = axis.binning();
        axis_data.bounds = axis.bounds();
        axis_data.label = axis.label();
        axis_data.label = axis.nbins();

        axis_data.edges.resize(2 * axis.nbins());
        std::copy(std::cbegin(axis.bin_edges()), std::cend(m.bin_edges()),
                  std::begin(axis_data.edges));

        return axis_data;
    }

    private:
    /// Retrieve @c mask_payload from mask_store element
    /*struct get_grid_payload {
        template <typename grid_coll_t, typename index_t>
        inline auto operator()(const grid_coll_t& grid_coll,
                               const index_t& index) const {
            return grid_writer<detector_t>::serialize(grid_coll[index]);
        }
    };*/
};

}  // namespace detray
