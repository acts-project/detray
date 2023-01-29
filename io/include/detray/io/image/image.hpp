/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracer/color.hpp"

// System include(s).
#include <iostream>
#include <vector>

namespace detray {

class image {

    public:
    /// Default constructor
    constexpr image() = default;

    /// Construct from image @param hight and @param width in pixels
    DETRAY_HOST_DEVICE
    image(const unsigned long hight, const unsigned long width,
          const color<float> c = {})
        : m_hight{hight}, m_width{width}, m_data{hight * width, c} {}

    /// @returns image hight in pixels
    DETRAY_HOST_DEVICE
    constexpr unsigned long hight() const { return m_hight; }

    /// @returns image width in pixels
    DETRAY_HOST_DEVICE
    constexpr unsigned long width() const { return m_width; }

    /// @returns number of pixels
    DETRAY_HOST_DEVICE
    constexpr unsigned long n_pixels() const { return m_hight * m_width; }

    /// @returns the pixel data
    DETRAY_HOST_DEVICE
    constexpr const std::vector<color<float>>& pixel_data() const {
        return m_data;
    }

    /// @returns the pixel data
    DETRAY_HOST_DEVICE
    constexpr std::vector<color<float>>& pixel_data() { return m_data; }

    private:
    /// Image size in pixels
    unsigned long m_hight{100u};
    unsigned long m_width{100u};

    /// Pixel data
    std::vector<color<float>> m_data{m_hight * m_width, color{}};
};

}  // namespace detray