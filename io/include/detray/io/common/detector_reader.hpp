/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/detector_components_io.hpp"
#include "detray/io/common/detail/type_traits.hpp"
#include "detray/io/json/json_reader.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/utils/consistency_checker.hpp"

// System include(s)
#include <filesystem>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace detray {

namespace io {

/// @brief config struct for detector reading.
struct detector_reader_config {
    /// Input files
    std::vector<std::string> m_files;
    /// Run detector consistency check after reading
    bool m_do_check{false};

    /// Getters
    /// @{
    const std::vector<std::string>& files() const { return m_files; }
    bool do_check() const { return m_do_check; }
    /// @}

    /// Setters
    /// @{
    detector_reader_config& add_file(const std::string file_name) {
        m_files.push_back(std::move(file_name));
        return *this;
    }
    detector_reader_config& do_check(const bool check) {
        m_do_check = check;
        return *this;
    }
    /// @}
};

}  // namespace io

namespace detail {

/// From the list of file that are given as part of the config @param cfg,
/// infer the readers that are needed by peeking into the file header
template <class detector_t>
auto assemble_reader(const io::detector_reader_config& cfg) {

    detail::detector_component_readers<detector_t> readers;

    // Read the name of the detector from file
    std::string det_name{""};

    for (const std::string& file_name : cfg.files()) {

        std::string extension{std::filesystem::path{file_name}.extension()};
        if (extension == ".json") {
            // Peek at the header to determine the kind of reader that is needed
            common_header_payload header = read_json_header(file_name);

            if (header.tag == "geometry") {
                det_name = header.detector;
                readers.template add<json_geometry_reader>(file_name);
            } else if (header.tag == "homogeneous_material") {
                readers.template add<json_homogeneous_material_reader>(
                    file_name);
            } else {
                throw std::invalid_argument("Unsupported file tag '" +
                                            header.tag +
                                            "' in input file: " + file_name);
            }
        } else {
            throw std::runtime_error("Unsupported file format '" + extension +
                                     "' for input file: " + file_name);
        }
    }

    return std::make_pair(std::move(readers), std::move(det_name));
}

}  // namespace detail

namespace io {

/// @brief Reader function for detray detectors.
///
/// @tparam detector_t the type of detector to be built
/// @tparam volume_builder_t the type of base volume builder to be used
///
/// @param resc the memory resource to be used for the detector container allocs
/// @param cfg the detector reader configuration
///
/// @returns a complete detector object + a map that contains the volume names
template <class detector_t,
          template <typename> class volume_builder_t = volume_builder>
auto read_detector(vecmem::memory_resource& resc,
                   const detector_reader_config& cfg) {
    // Map the volue names to their indices
    typename detector_t::name_map names{};

    detector_builder<typename detector_t::metadata, volume_builder_t>
        det_builder;

    // Find all required the readers
    auto [reader, det_name] = detray::detail::assemble_reader<detector_t>(cfg);

    names.emplace(0u, std::move(det_name));
    reader.read(det_builder, names);

    // Build and return the detector
    auto det = det_builder.build(resc);

    if (cfg.do_check()) {
        detray::detail::check_consistency(det);
    }

    return std::make_pair(std::move(det), std::move(names));
}

}  // namespace io

}  // namespace detray
