/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <array>
#include <string>
#include <vector>

namespace detray::io {

/// @brief config struct for detector reading.
struct detector_reader_config {
    /// Input files
    std::vector<std::string> m_files;
    /// Run detector consistency check after reading
    bool m_do_check{true};

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

}  // namespace detray::io
