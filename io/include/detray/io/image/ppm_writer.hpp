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
template <typename color_depth = uint8_t>
class ppm_writer final : public image_writer<color_depth> {

    public:
    /// All writers must define a file name
    ppm_writer() = delete;

    /// File gets created with a @param name and @param extension
    ppm_writer(const std::string &name = "")
        : image_writer<color_depth>(name, "ppm") {}

    /// Writes the image to disk as a ppm file
    void write(const raw_image<color_depth> &im) override {
        // ppm file header
        this->m_file << "P3\n" << im.width() << " " << im.hight() << "\n255\n";

        // Image data
        for (const auto &px : im.pixel_data()) {
            this->m_file << static_cast<uint>(px[0]) << " "
                         << static_cast<uint>(px[1]) << " "
                         << static_cast<uint>(px[2]) << std::endl;
        }
    };
};

}  // namespace detray::io
