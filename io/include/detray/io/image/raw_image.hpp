/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/pixel.hpp"

// System include(s).
#include <iostream>
#include <vector>

namespace detray::io {

template <typename depth = uint8_t>
class raw_image {

    using color_t = texture::color<depth>;

    public:
    /// Default constructor
    constexpr raw_image() = default;

    /// Construct from image @param height and @param width in pixels
    DETRAY_HOST_DEVICE
    raw_image(const unsigned int width, const unsigned int height,
              const color_t c = {})
        : m_width{width}, m_height{height}, m_data{height * width, c} {}

    /// @returns image width in pixels
    DETRAY_HOST_DEVICE
    constexpr unsigned long width() const { return m_width; }

    /// @returns image height in pixels
    DETRAY_HOST_DEVICE
    constexpr unsigned long height() const { return m_height; }

    /// @returns number of pixels
    DETRAY_HOST_DEVICE
    constexpr unsigned long n_pixels() const { return m_height * m_width; }

    /// @returns the pixel data
    DETRAY_HOST_DEVICE
    constexpr const std::vector<color_t>& pixel_data() const { return m_data; }

    /// @returns the pixel data
    DETRAY_HOST_DEVICE
    constexpr std::vector<color_t>& pixel_data() { return m_data; }

    /// Set a particular pixel in the image
    DETRAY_HOST_DEVICE
    constexpr void set_pixel(uint x, uint y, color_t c) {
        std::size_t px_idx{x + m_width * y};
        m_data.at(px_idx) = c;
    }

    /// Set a particular pixel in the image to @param px
    DETRAY_HOST_DEVICE
    constexpr void set_pixel(texture::pixel<unsigned int, depth>& px) {
        set_pixel(px[0], px[1], px.color());
    }

    private:
    /// Image size in pixels
    unsigned long m_width{100u};
    unsigned long m_height{100u};

    /// Pixel data
    std::vector<color_t> m_data{m_height * m_width, color_t{}};
};

}  // namespace detray::io
