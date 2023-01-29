/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/io/image/image.hpp"

// System include(s).
#include <fstream>
#include <iostream>
#include <string>

namespace detray {

/// @brief Abstract base class for image writers
class image_writer {

    public:
    /// All writers must define a file name
    image_writer() = delete;

    /// File gets created with a @param name and @param extension
    image_writer(const std::string &name = "",
                 const std::string &extension = "img") {
        // default name
        std::string file_name{
            name.empty() ? "graphic_" + std::to_string(n_files) : name};

        // Count the new file
        ++n_files;

        // Open file
        m_file.open(file_name + "." + extension);
    }

    /// Destructor closes the file
    virtual ~image_writer() { m_file.close(); }

    /// Writes an image to disk
    virtual void write(const image &im) = 0;

    protected:
    /// Output file handle
    std::ofstream m_file;

    private:
    /// How many writers/files have been created?
    static unsigned int n_files;
};

unsigned int image_writer::n_files{0u};

}  // namespace detray
