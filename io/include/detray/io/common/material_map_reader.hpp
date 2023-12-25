/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/type_traits.hpp"
#include "detray/io/common/homogeneous_material_reader.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/tools/material_map_builder.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <memory>
#include <stdexcept>
#include <vector>

namespace detray {

/// @brief Abstract base class for material map readers
template <class detector_t,
          typename DIM = std::integral_constant<std::size_t, 2u>>
class material_map_reader : public reader_interface<detector_t> {

    using scalar_t = typename detector_t::scalar_type;
    using mat_id = typename detector_t::materials::id;
    using material_reader_t = homogeneous_material_reader<detector_t>;

    public:
    static constexpr std::size_t dim{DIM()};

    using base_type = reader_interface<detector_t>;
    using bin_index_type = axis::multi_bin<dim>;

    /// Same constructors for this class as for base_type
    using base_type::base_type;

    private:
    using mat_factory_t = material_map_factory<detector_t, bin_index_type>;
    using mat_data_t = typename mat_factory_t::data_type;
    /// Gets compile-time mask information
    template <io::detail::material_type mat_type>
    using map_info = detail::mat_map_info<mat_type, detector_t>;

    protected:
    /// Tag the reader as "material_maps"
    inline static const std::string tag = "material_maps";

    /// Deserialize the material grids @param grids_data from their IO
    /// payload
    static void deserialize(
        detector_builder<typename detector_t::metadata, volume_builder>
            &det_builder,
        typename detector_t::name_map &,
        detector_grids_payload<material_slab_payload, io::detail::material_type>
            &&grids_data) {

        // Deserialize the material volume by volume
        for (const auto &[vol_idx, mat_grids] : grids_data.grids) {
            // Decorate the current volume builder with material maps
            auto vm_builder =
                det_builder
                    .template decorate<material_map_builder<detector_t, dim>>(
                        static_cast<dindex>(vol_idx));

            // Add the material data to the factory
            auto mat_factory = std::make_shared<
                material_map_factory<detector_t, bin_index_type>>();

            // Deserialize the material grid of each surface
            for (const auto &grid_data : mat_grids) {

                mat_id map_id = deserialize<io::detail::material_type::n_mats>(
                    grid_data.grid_link.type);

                // Get the number of bins per axis
                std::vector<std::size_t> n_bins{};
                for (const auto &axis_data : grid_data.axes) {
                    n_bins.push_back(axis_data.bins);
                }

                // Get the local bin indices and the material parametrization
                std::vector<bin_index_type> loc_bins{};
                mat_data_t mat_data{
                    base_type::deserialize(grid_data.owner_link)};
                for (const auto &bin_data : grid_data.bins) {

                    assert(dim == bin_data.loc_index.size() &&
                           "Dimension of local bin indices in input file does "
                           "not match material grid dimension");

                    // The local bin indices for the bin to be filled
                    bin_index_type mbin;
                    for (const auto &[i, bin_idx] :
                         detray::views::enumerate(bin_data.loc_index)) {
                        mbin[i] = bin_idx;
                    }
                    loc_bins.push_back(std::move(mbin));

                    // For now assume surfaces ids as the only grid input
                    for (const auto &slab_data : bin_data.content) {
                        mat_data.append(
                            material_reader_t::deserialize(slab_data));
                    }
                }

                mat_factory->add_material(map_id, std::move(mat_data),
                                          std::move(n_bins),
                                          std::move(loc_bins));
            }

            // Add the material maps to the volume
            vm_builder->add_surfaces(mat_factory);
        }
    }

    private:
    /// Get the detector material id from the payload material type id
    template <io::detail::material_type I>
    static mat_id deserialize(io::detail::material_type type_id) {

        // Material id of map data found
        if (type_id == I) {
            // Get the corresponding material id for this detector
            return map_info<I>::value;
        }
        // Test next material type id
        constexpr int current_id{static_cast<int>(I)};
        if constexpr (current_id > 0) {
            return deserialize<static_cast<io::detail::material_type>(
                current_id - 1)>(type_id);
        } else {
            return mat_id::e_none;
        }
    }
};

}  // namespace detray
