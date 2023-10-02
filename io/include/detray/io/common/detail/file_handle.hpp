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
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace detray::io {

enum class format { json = 0u };

namespace detail {

/// @brief Wrapper around a file stream
///
/// Performs some checks and offers some convenience functionality:
/// - In out mode w.o. file replacement: Checks if file exists and automatically
///   finds a new filename in case it does.
/// - In in mode: Checks for empty file names and whether the file exists
/// - Enforces a limit on the number of files that can be written/opened
/// - Checks whether a file was opened correctly
/// - Closes the stream when the handle goes out of scope and checks whether
///   anything went wrong during the IO operations
///
/// @note Not thread safe.
/// @note Can throw exceptions during construction.
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
        // File name
        std::string file_stem{name};

        // Pure output mode without replacement of file: Check if name is taken
        // and modify it if necessary
        if (mode == std::ios_base::out or
            (mode == (std::ios_base::out | std::ios_base::binary))) {
            // Default name for output
            file_stem =
                name.empty() ? "detray_" + std::to_string(n_files) : name;

            std::filesystem::path file_path{file_stem + extension};
            std::size_t n_trials{2u};
            while (std::filesystem::exists(file_path)) {
                file_stem = get_alternate_file_stem(file_stem, n_trials);
                file_path = std::filesystem::path{file_stem + extension};
                ++n_trials;

                // The maximum here is arbitrary
                if (n_trials >= 10000u) {
                    throw std::runtime_error(
                        "Too many versions of file exist: " + file_stem +
                        extension);
                }
            }
            // Pure input mode: Check if file name makes sense and file exists
        } else if ((mode == std::ios_base::in) or
                   (mode == (std::ios_base::in | std::ios_base::binary))) {
            if (file_stem.empty()) {
                throw std::invalid_argument("File name empty");
            }

            std::filesystem::path file_path{file_stem + extension};
            if (not std::filesystem::exists(file_path)) {
                throw std::invalid_argument(
                    "Could not open file: File does not exist: " + file_stem +
                    extension);
            }
        } else if (file_stem.empty()) {
            std::cout << "WARNING: Empty file name" << std::endl;
        }

        // Count the new file
        const std::string file_name{file_stem + extension};
        ++n_files;
        ++n_open_files;
        if (n_files >= std::numeric_limits<std::uint_least16_t>::max()) {
            throw std::runtime_error(
                "Could not open file: Too many files written: " + file_name);
        } else if (n_open_files >= 1000u) {
            throw std::runtime_error(
                "Could not open file: Too many files currently open: " +
                file_name);
        }

        // Open file
        m_stream.open(file_name, mode);

        if (!m_stream.is_open()) {
            throw std::runtime_error("Could not open file: " + file_name);
        }
    }

    /// Destructor closes the file
    ~file_handle() {
        if (m_stream.bad()) {
            std::cout << "ERROR: Could not read from/write to file";
        }
        m_stream.close();
        --n_open_files;
    }

    /// @returns the output stream
    std::fstream& operator*() { return m_stream; }

    private:
    /// @returns alternate file stem upon collision
    std::string get_alternate_file_stem(std::string& stem,
                                        const std::size_t n) {
        const std::string delim{"_"};

        // File stem already comes with a number, simply update it
        if (n > 2u) {
            std::size_t pos{stem.rfind(delim)};
            if (pos == std::string::npos or (pos + 1 == stem.size())) {
                throw std::runtime_error("Malformed file name");
            }

            // Leave the delimeter where it is and replace the file number
            ++pos;
            stem.replace(pos, stem.size() - pos, std::to_string(n));

            return stem;
        }

        return stem + delim + std::to_string(n);
    }

    /// Output file handle
    std::fstream m_stream;

    /// How many files have been created? Maximum: 65'536
    inline static std::size_t n_files{0u};
    inline static std::size_t n_open_files{0u};
};

}  // namespace detail

}  // namespace detray::io
