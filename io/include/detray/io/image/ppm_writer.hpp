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
template <typename color_depth = std::uint8_t>
class ppm_writer final : public image_writer<color_depth> {

    public:
    /// File gets created with a @param name and @param extension
    ppm_writer() : image_writer<color_depth>(".ppm") {}

    /// Writes the image to disk as a ppm file
    void write(const raw_image<color_depth> &im,
               const std::string &file_stem) override {
        // Create a new file
        io::detail::file_handle file{file_stem, this->m_file_extension, (std::ios_base::out | std::ios_base::trunc)};

        // ppm file header
        *file << "P3\n" << im.width() << " " << im.height() << "\n255\n";

        // Image data
        for (const auto &px : im.pixel_data()) {
            *file << static_cast<uint>(px[0]) << " " << static_cast<uint>(px[1])
                  << " " << static_cast<uint>(px[2]) << std::endl;
        }
    };
};

}  // namespace detray::io
