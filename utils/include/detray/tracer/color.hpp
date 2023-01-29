/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/definitions/math.hpp"

// System include(s)
#include <array>
#include <cmath>
#include <iostream>

namespace detray {

/// @brief holds rgb and alpha values for color shading
///
/// @tparam data_t how to store the rgb data: single color or soa
template<typename data_t = float>
struct color {

    /// Default constructor
    constexpr color() = default;

    /// Construct from colors @param r (red), @param g (green), @param b (blue)
    /// and @param alpha values
    DETRAY_HOST_DEVICE
    constexpr color(const data_t r, const data_t g, const data_t b, const data_t alpha) 
        : m_data{r, g, b, alpha}{}

    /// Subscript operator @returns a color data point - const
    DETRAY_HOST_DEVICE
    constexpr data_t operator[](const std::size_t i) const {
        return m_data[i];
    }

    /// Subscript operator @returns a color data point - non-const
    DETRAY_HOST_DEVICE
    constexpr data_t operator[](const std::size_t i) {
        return m_data[i];
    }

    /// Print the color data to stdout
    DETRAY_HOST
    template<typename> 
    friend std::ostream & operator<<(std::ostream &os, const color& c);

    std::array<data_t, 4> m_data{};
};

template<typename data_t>
std::ostream & operator<<(std::ostream &os, const color<data_t>& c) {
    return os << "rgba: (" << c[0] << ", " << c[1] << ", " << c[2] << ", " << c[3] << ")";
}

}  // namespace detray
