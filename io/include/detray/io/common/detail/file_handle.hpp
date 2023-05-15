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
#include <filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <stdexcept>
#include <string>

namespace detray::io {

enum class format { json = 0u };

namespace detail {

/// Wrapper around a file stream
class file_handle final {

    public:
    /// All writers must define a file name
    file_handle() = delete;

    /// File gets created with a @param name and @param extension
    file_handle(const std::string& file_name,
                std::ios_base::openmode mode = std::ios_base::in |
                                               std::ios_base::out)
        : file_handle(std::filesystem::path{file_name}.parent_path() /
                          std::filesystem::path{file_name}.stem(),
                      std::filesystem::path{file_name}.extension(), mode) {}

    /// File gets created with a @param name and @param extension
    file_handle(const std::string& name, const std::string& extension,
                std::ios_base::openmode mode = std::ios_base::in |
                                               std::ios_base::out) {
        // Default name
        std::string file_stem{name.empty() ? "detray_" + std::to_string(n_files)
                                           : name};

        // Check if name is taken and modify it if necessary
        if (mode == std::ios_base::out) {
            std::filesystem::path file_path{file_stem + extension};
            std::size_t n_trials{1u};
            while (std::filesystem::exists(file_path)) {
                file_stem = get_alternate_file_stem(file_stem, n_trials);
                file_path = std::filesystem::path{file_stem + extension};
                ++n_trials;
                // The maximum here is arbitrary
                if (n_trials >= 10000u) {
                    throw std::runtime_error("Too many versions of file exist");
                }
            }
        }

        // Count the new file
        ++n_files;
        ++n_open_files;
        assert(n_files < std::numeric_limits<std::uint_least16_t>::max());
        assert(n_open_files < 1000u);

        // Open file
        m_stream.open(file_stem + extension, mode);

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
    /// @returns alternate file stem upon collision
    std::string get_alternate_file_stem(std::string& stem,
                                        const std::size_t n) {
        // File stem already come with a number, simply update it
        if (n > 1u) {
            stem.replace(stem.size() - 1u, 1, std::to_string(n));
            return stem;
        }
        return stem + "_" + std::to_string(n);
    }

    /// Output file handle
    std::fstream m_stream;

    /// How many files have been created? Maximum: 65'536
    static std::size_t n_files;
    static std::size_t n_open_files;
};

std::size_t file_handle::n_files{0u};
std::size_t file_handle::n_open_files{0u};

}  // namespace detail

}  // namespace detray::io
