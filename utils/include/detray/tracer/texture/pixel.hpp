/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracer/texture/color.hpp"

// System include(s)
#include <iostream>

namespace detray::texture {

namespace detail {

/// @brief holds rgb and alpha values for color shading
///
/// @tparam data_t how to store the rgb data: single color or soa
template <uint D, typename data_t = uint, typename depth = uint8_t>
struct pixelD {

    using color_t = texture::color<depth>;

    static constexpr uint Dim{D};

    /// Default constructor
    constexpr pixelD() = default;

    /// Construct from colors @param r (red), @param g (green), @param b (blue)
    /// and @param alpha values
    DETRAY_HOST_DEVICE
    constexpr pixelD(const std::array<data_t, D>& coord, const color_t& c)
        : m_color{c}, m_coord{coord} {}

    /// Equality operator: Only considers exact match
    DETRAY_HOST_DEVICE
    constexpr data_t operator==(const pixelD& other) {
        return (m_coord == other.m_coord) and (m_color == other.m_color);
    }

    /// Subscript operator @returns a color data point - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator[](const std::size_t i) const {
        return m_coord[i];
    }

    /// Subscript operator @returns a color data point - non-const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator[](const std::size_t i) {
        return m_coord[i];
    }

    /// @returns the color of the pixel - non-const
    DETRAY_HOST_DEVICE
    constexpr const color_t& color() const { return m_color; }

    /// Set the color of the pixel to @param c
    DETRAY_HOST_DEVICE
    constexpr void set_color(const color_t& c) { m_color = c; }

    /// Mixes the by averaging their positions and mixing their colors
    DETRAY_HOST_DEVICE
    template <uint, typename, typename>
    friend constexpr pixelD operator+(const pixelD&, const pixelD&);

    /// Print the pixel data to stdout
    DETRAY_HOST
    template <uint, typename, typename>
    friend std::ostream& operator<<(std::ostream&, const pixelD&);

    color_t m_color{};
    std::array<data_t, D> m_coord{};
};

template <uint D, typename data_t, typename depth>
std::ostream& operator<<(std::ostream& os, const pixelD<D, data_t, depth>& px) {
    return os << "pix: " << static_cast<uint>(px[0]) << ", "
              << static_cast<uint>(px[1]) << ", " << px.color();
}

template <uint D, typename data_t, typename depth>
constexpr pixelD<D, data_t, depth> operator+(
    const pixelD<D, data_t, depth>& left,
    const pixelD<D, data_t, depth>& right) {
    pixelD<D, data_t, depth> new_px;
    new_px.set_color(left.color() + right.color());

    for (uint i{0u}; i < pixelD<D, data_t, depth>::Dim; ++i) {
        new_px[i] = (left[i] + right[i]) / 2u;
    }

    return new_px;
}

}  // namespace detail

template <typename data_t = uint, typename depth = uint8_t>
using pixel = detail::pixelD<2, data_t, depth>;

template <typename data_t = uint, typename depth = uint8_t>
using voxel = detail::pixelD<3, data_t, depth>;

}  // namespace detray::texture
