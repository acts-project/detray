/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/io/frontend/detail/detector_components_reader.hpp"
#include "detray/io/frontend/detector_reader_config.hpp"
#include "detray/io/frontend/impl/json_readers.hpp"
#include "detray/utils/consistency_checker.hpp"

// System include(s)
#include <filesystem>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace detray::io {

/// @brief Read detector components from a list of files
///
/// @param file_names list of files to be read
/// @param det_builder detector builder to be filled
/// @param name_map detector and volume name map
template <class detector_t, std::size_t CAP = 0u, std::size_t DIM = 2u>
void read_components_from_file(const std::vector<std::string>& file_names,
                               detector_builder<typename detector_t::metadata,
                                                volume_builder>& det_builder,
                               typename detector_t::name_map& name_map) {
    // Hold all required readers (one for every component)
    detail::detector_components_reader<detector_t> readers;

    // Register the readers for the files in json format
    detail::add_json_readers<CAP, DIM>(readers, file_names);

    // Make sure that all files will be read
    if (readers.size() != file_names.size()) {
        std::cout << "WARNING: Not all files were registered to a reader. "
                  << "Registered files: ";

        for (const auto& [file_name, reader] : readers.readers_map()) {
            std::cout << "\t" << file_name << std::endl;
        }
    }

    // Read the data into the detector builder
    readers.read(det_builder, name_map);
}

/// @brief Reader function for detray detectors.
///
/// @tparam detector_t the type of detector to be built
/// @tparam CAP surface grid bin capacity. If CAP is 0, the grid reader builds
///             a grid type with dynamic bin capacity
/// @tparam DIM dimension of the surface grids, usually 2D
/// @tparam volume_builder_t the type of base volume builder to be used
///
/// @param resc the memory resource to be used for the detector container allocs
/// @param cfg the detector reader configuration
///
/// @returns a complete detector object + a map that contains the volume names
template <class detector_t, std::size_t CAP = 0u, std::size_t DIM = 2u,
          template <typename> class volume_builder_t = volume_builder>
auto read_detector(vecmem::memory_resource& resc,
                   const detector_reader_config& cfg) noexcept(false) {

    // Map the volume names to their indices
    typename detector_t::name_map names{};

    detector_builder<typename detector_t::metadata, volume_builder_t>
        det_builder;

    // Register readers for the respective detector component and file format
    // and read the data into the detector_builder
    read_components_from_file<detector_t, CAP, DIM>(cfg.files(), det_builder,
                                                    names);

    // Build and return the detector
    auto det = det_builder.build(resc);

    if (cfg.do_check()) {
        // This will throw an exception in case of inconsistencies
        detray::detail::check_consistency(det, cfg.verbose_check(), names);
        std::cout << "Detector check: OK" << std::endl;
    }

    return std::make_pair(std::move(det), std::move(names));
}

}  // namespace detray::io
