/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/bfield_reader.hpp"
#include "detray/io/common/detail/detector_components_io.hpp"
#include "detray/io/common/detail/type_traits.hpp"
#include "detray/io/json/json_reader.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/utils/consistency_checker.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/vector.hpp>

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
    bool m_do_check{true};
    /// Field vector for an homogenoues b-field
    std::array<scalar, 3> m_bfield_vec{};

    /// Getters
    /// @{
    const std::vector<std::string>& files() const { return m_files; }
    bool do_check() const { return m_do_check; }
    template <typename const_bfield_bknd_t>
    auto bfield() const {
        return covfie::field<const_bfield_bknd_t>(covfie::make_parameter_pack(
            typename const_bfield_bknd_t::configuration_t{
                m_bfield_vec[0], m_bfield_vec[1], m_bfield_vec[2]}));
    }
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
    detector_reader_config& bfield_vec(const scalar x, const scalar y,
                                       const scalar z) {
        m_bfield_vec = {x, y, z};
        return *this;
    }
    /// @}
};

}  // namespace io

namespace detail {

/// From the list of files that are given as part of the config @param cfg,
/// infer the readers that are needed by peeking into the file headers
///
/// @tparam detector_t type of the detector instance: Must match the data that
///                    is read from file!
/// @tparam CAP surface grid bin capacity (@TODO make runtime)
/// @tparam DIM dimension of the surface grids, usually 2D
template <class detector_t, std::size_t CAP, std::size_t DIM>
auto assemble_reader(const io::detector_reader_config& cfg) noexcept(false) {

    detail::detector_component_readers<detector_t> readers;

    // Read the name of the detector from file
    std::string det_name{""};

    for (const std::string& file_name : cfg.files()) {

        if (file_name.empty()) {
            std::cout << "WARNING: Empty file name. Component will not be built"
                      << std::endl;
            continue;
        }

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

            } else if (header.tag == "surface_grids") {
                using surface_t = typename detector_t::surface_type;
                readers.template add<json_grid_reader, surface_t,
                                     std::integral_constant<std::size_t, CAP>,
                                     std::integral_constant<std::size_t, DIM>>(
                    file_name);
            } else {
                throw std::invalid_argument("Unsupported file tag '" +
                                            header.tag +
                                            "' in input file: " + file_name);
            }

            // This is the file type covfie uses
        } else if (extension == ".cvf" and check_covfie_file(file_name)) {

            if constexpr (not std::is_same_v<
                              typename detector_t::bfield_type::backend_t,
                              bfield::const_bknd_t>) {
                readers.template add<covfie_reader>(file_name);
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
/// @tparam CAP surface grid bin capacity (@TODO make runtime)
/// @tparam DIM dimension of the surface grids, usually 2D
/// @tparam volume_builder_t the type of base volume builder to be used
///
/// @param resc the memory resource to be used for the detector container allocs
/// @param cfg the detector reader configuration
///
/// @returns a complete detector object + a map that contains the volume names
template <class detector_t, std::size_t CAP = 9u, std::size_t DIM = 2u,
          template <typename> class volume_builder_t = volume_builder>
auto read_detector(vecmem::memory_resource& resc,
                   const detector_reader_config& cfg) noexcept(false) {

    using bfield_bknd_t = typename detector_t::bfield_type::backend_t;

    // Map the volume names to their indices
    typename detector_t::name_map names{};

    detector_builder<typename detector_t::metadata, bfield_bknd_t,
                     volume_builder_t>
        det_builder;

    // Find all required readers
    auto [reader, det_name] =
        detray::detail::assemble_reader<detector_t, CAP, DIM>(cfg);

    // Add a constant b-field (no bfield reader is added in this case)
    if constexpr (std::is_same_v<bfield_bknd_t, bfield::const_bknd_t>) {
        det_builder.set_bfield(cfg.template bfield<bfield_bknd_t>());
    }

    // Read the data
    names.emplace(0u, std::move(det_name));
    reader.read(det_builder, names);

    // Build and return the detector
    auto det = det_builder.build(resc);

    if (cfg.do_check()) {
        // This will throw an exception in case of inconsistencies
        detray::detail::check_consistency(det);
        std::cout << "Detector check: OK" << std::endl;
    }

    return std::make_pair(std::move(det), std::move(names));
}

}  // namespace io

}  // namespace detray
