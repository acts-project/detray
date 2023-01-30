/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/io/image/raw_image.hpp"

// System include(s)
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

namespace detray::io {

namespace detail {

/// Wrapper around a file stream for an image file
class file_handle final {

    public:
    /// All writers must define a file name
    file_handle() = delete;

    /// File gets created with a @param name and @param extension
    file_handle(const std::string &name = "",
                const std::string &extension = "img") {
        // Default name
        std::string file_name{
            name.empty() ? "graphic_" + std::to_string(n_files) : name};

        // Count the new file
        ++n_files;
        assert(n_files < std::numeric_limits<uint>::max());

        // Open file
        m_out_stream.open(file_name + "." + extension);
    }

    /// Destructor closes the file
    ~file_handle() { m_out_stream.close(); }

    /// @returns the output stream
    std::ofstream &operator*() { return m_out_stream; }

    protected:
    /// Output file handle
    std::ofstream m_out_stream;

    private:
    /// How many writers/files have been created?
    static unsigned int n_files;
};

unsigned int file_handle::n_files{0u};

}  // namespace detail

/// @brief Abstract base class for image writers
template <typename color_depth = uint8_t>
class image_writer {

    public:
    /// All writers must define a file name
    image_writer() = delete;

    /// File gets created with a fixed @param extension
    image_writer(const std::string &ext = "img") : m_file_extension{ext} {}

    /// Default destructor
    virtual ~image_writer() {}

    /// Writes an image to file with a given name
    virtual void write(const raw_image<color_depth> &, const std::string &) = 0;

    protected:
    std::string m_file_extension;
};

}  // namespace detray::io
