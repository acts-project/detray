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
#include "detray/io/common/detail/utils.hpp"
#include "detray/io/json/json_writer.hpp"

// System include(s)
#include <filesystem>
#include <ios>

namespace detray {

namespace io {

/// @brief config struct for detector writing.
struct detector_writer_config {
    /// The path to the output files
    std::string m_path{"./"};
    /// The output file format
    detray::io::format m_format = detray::io::format::json;
    /// Replace files in case they already exists
    bool m_replace = false;
    /// Compactify json output, if not json format this flag does nothing
    bool m_compact_io = false;
    /// Whether to write the material to file
    bool m_write_material = true;
    /// Whether to write the accelerator grids to file
    bool m_write_grids = true;

    /// Getters
    /// @{
    const std::string& path() const { return m_path; }
    detray::io::format format() const { return m_format; }
    bool replace_files() const { return m_replace; }
    bool compactify_json() const { return m_compact_io; }
    bool write_material() const { return m_write_material; }
    bool write_grids() const { return m_write_grids; }
    /// @}

    /// Setters
    /// @{
    detector_writer_config& path(std::string p) {
        m_path = std::move(p);
        return *this;
    }
    detector_writer_config& format(detray::io::format f) {
        m_format = f;
        return *this;
    }
    detector_writer_config& replace_files(bool flag) {
        m_replace = flag;
        return *this;
    }
    detector_writer_config& compactify_json(bool flag) {
        m_compact_io = flag;
        return *this;
    }
    detector_writer_config& write_material(bool flag) {
        m_write_material = flag;
        return *this;
    }
    detector_writer_config& write_grids(bool flag) {
        m_write_grids = flag;
        return *this;
    }
    /// @}
};

}  // namespace io

namespace detail {

/// From the detector type @tparam detector_t, infer the writers that are needed
template <class detector_t>
detail::detector_component_writers<detector_t> assemble_writer(
    const io::detector_writer_config& cfg) {

    detail::detector_component_writers<detector_t> writers;

    if (cfg.format() == io::format::json) {
        // Always needed
        writers.template add<json_geometry_writer>();

        // Find other writers, depending on the detector type
        if (cfg.write_material()) {
            // Simple material
            if constexpr (detail::is_homogeneous_material_v<detector_t>) {
                writers.template add<json_homogeneous_material_writer>();
            }
            // Material maps
            // ...
        }
        if (cfg.write_grids()) {
            writers.template add<json_grid_writer>();
        }
    }

    return writers;
}

}  // namespace detail

namespace io {

/// @brief Writer function for detray detectors.
///
/// Based on both the given config/file format and the detector type,
/// the correct writers are being assembled and called.
/// Write @param det to file in the format @param format,
/// using @param names to name the components
template <class detector_t>
void write_detector(detector_t& det, const typename detector_t::name_map& names,
                    const detector_writer_config& cfg) {
    // How to open the file
    const std::ios_base::openmode out_mode{std::ios_base::out |
                                           std::ios_base::binary};
    auto mode =
        cfg.replace_files() ? (out_mode | std::ios_base::trunc) : out_mode;

    const auto file_path = detray::detail::create_path(cfg.path());

    // Get the writer
    auto writer = detray::detail::assemble_writer<detector_t>(cfg);

    writer.write(det, names, mode, file_path);
}

}  // namespace io

}  // namespace detray
