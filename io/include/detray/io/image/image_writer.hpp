/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/tracer/renderer/raw_image.hpp"
#include "detray/io/common/detail/file_handle.hpp"

// System include(s)
#include <string>

namespace detray::io {

/// @brief Abstract base class for image writers
template <typename color_depth = std::uint8_t>
class image_writer {

    public:
    /// All writers must define a file name
    image_writer() = delete;

    /// File gets created with a fixed @param extension
    image_writer(const std::string &ext = ".img") : m_file_extension{ext} {}

    /// Default destructor
    virtual ~image_writer() {}

    /// Writes an image to file with a given name
    virtual void write(const raw_image<color_depth> &, const std::string &) = 0;

    protected:
    std::string m_file_extension;
};

}  // namespace detray::io
