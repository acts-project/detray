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
#include <ratio>
#include <vector>

namespace detray {

template <typename depth = std::uint8_t, typename ratio_t = std::ratio<16, 9>>
class raw_image {

    public:
    using color_t = texture::color<depth>;
    using aspect_ratio = ratio_t;

    /// Default constructor
    constexpr raw_image() = default;

    /// Construct from image @param height and @param width in pixels
    DETRAY_HOST_DEVICE
    raw_image(const unsigned int height, const color_t c = {})
        : m_width{static_cast<std::size_t>(
              height * static_cast<double>(aspect_ratio::num) /
              static_cast<double>(aspect_ratio::den))},
          m_height{height},
          m_data{m_height * m_width, c} {}

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
    constexpr void set_pixel(std::uint64_t x, std::uint64_t y, color_t c) {
        std::size_t px_idx{x + m_width * y};
        m_data.at(px_idx) = c;
    }

    /// Set a particular pixel in the image to @param px
    DETRAY_HOST_DEVICE
    constexpr void set_pixel(texture::pixel<std::uint64_t, depth>& px) {
        set_pixel(px[0], px[1], px.color());
    }

    private:
    /// Image size in pixels
    std::size_t m_width{160u};
    std::size_t m_height{90u};

    /// Pixel data
    std::vector<color_t> m_data{m_height * m_width, color_t{}};
};

}  // namespace detray
