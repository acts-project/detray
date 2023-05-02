/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <cassert>
#include <cstdint>
#include <fstream>
#include <ios>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace detray::io::detail {

/// Wrapper around a file stream
class file_handle final {

    public:
    /// All writers must define a file name
    file_handle() = delete;

    /// File gets created with a @param name and @param extension
    file_handle(const std::string& name = "",
                const std::string& extension = "txt",
                std::ios_base::openmode mode = std::ios_base::in |
                                               std::ios_base::out) {
        // Default name
        std::string file_name{name.empty() ? "detray_" + std::to_string(n_files)
                                           : name};

        // Count the new file
        ++n_files;
        ++n_open_files;
        assert(n_files < std::numeric_limits<std::uint_least16_t>::max());
        assert(n_open_files < 1000u);

        // Open file
        m_stream.open(file_name + "." + extension, mode);

        if (!m_stream.is_open()) {
            throw std::runtime_error("Could not open file");
        }
    }

    /// Destructor closes the file
    ~file_handle() {
        m_stream.close();
        --n_open_files;
    }

    /// @returns the output stream
    std::fstream& operator*() { return m_stream; }

    private:
    /// Output file handle
    std::fstream m_stream;

    /// How many files have been created? Maximum: 65'536
    static std::size_t n_files;
    static std::size_t n_open_files;
};

std::size_t file_handle::n_files{0u};
std::size_t file_handle::n_open_files{0u};

}  // namespace detray::io::detail
