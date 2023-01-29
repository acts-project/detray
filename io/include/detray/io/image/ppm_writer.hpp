/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/io/image/image_writer.hpp"

namespace detray::io {

/// @brief This class is a file writer for the ppm format
class ppm_writer final : public image_writer {

    public:
    /// All writers must define a file name
    ppm_writer() = delete;

    /// File gets created with a @param name and @param extension
    ppm_writer(const std::string &name = "") : image_writer(name, "ppm") {}

    /// Writes the image to disk as a ppm file
    void write(const raw_image &im) override {
        // ppm file header
        m_file << "P3\n" << im.width() << " " << im.hight() << "\n255\n";

        // Image data
        for (const auto &px : im.pixel_data()) {
            m_file << static_cast<uint>(px[0]) << " "
                   << static_cast<uint>(px[1]) << " "
                   << static_cast<uint>(px[2]) << std::endl;
        }
    };
};

}  // namespace detray::io
