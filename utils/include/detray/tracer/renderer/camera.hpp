/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/tracer/renderer/raw_image.hpp"

// System include(s)
#include <ratio>

namespace detray {

template <typename scalar_t, typename aspect_ratio = std::ratio<16, 9>>
class camera {
    public:
    DETRAY_HOST_DEVICE
    camera(const scalar_t viewport_height = 2.0f,
           const point3D origin = {0.0f, 0.0f, 0.0f},
           const scalar_t focal_length = 1.0f)
        : m_origin(origin) {
        constexpr scalar_t a{static_cast<scalar_t>(aspect_ratio::num) /
                             static_cast<scalar_t>(aspect_ratio::den)};
        const scalar_t viewport_width{a * viewport_height};

        m_horizontal = point3D{viewport_width, 0.f, 0.f};
        m_vertical = point3D{0.f, viewport_height, 0.f};
        point3D av = (m_horizontal + m_vertical);
        m_lower_left_corner = m_origin - 0.5f * av -
                              vector3D{0.f, 0.f, focal_length};
    }

    DETRAY_HOST_DEVICE
    template <typename color_depth>
    constexpr ray<transform3D> get_ray(
        const scalar_t x, const scalar_t y,
        const raw_image<color_depth> &image) const {

        // percentage of pixel position of the width/height of the image
        const float u{x / static_cast<scalar_t>(image.width() - 1u)};
        const float v{y / static_cast<scalar_t>(image.height() - 1u)};

        return {
            m_origin,
            m_lower_left_corner + u * m_horizontal + v * m_vertical - m_origin,
            0.f, 0.f};
    }

    private:
    point3D m_origin;
    point3D m_lower_left_corner;
    vector3D m_horizontal;
    vector3D m_vertical;
};

}  // namespace detray
