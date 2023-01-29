/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/pixel.hpp"

// System include(s)
#include <iostream>

namespace detray {

constexpr texture::color<> red{139u, 0u, 0u, 0u};
constexpr texture::color<> blue{0u, 0u, 139u, 0u};
constexpr texture::color<> purple{red + blue};

/// Calculates the color of a pixel. Starting point of the shader pipeline
template <typename pixel_coord_t = uint, typename color_depth = uint8_t>
struct colorizer : detray::actor {

    struct state {

        state(const pixel_coord_t x, const pixel_coord_t y) : m_pixel{{x, y}} {}
        /// The resulting pixel
        texture::pixel<pixel_coord_t, color_depth> m_pixel;
    };

    /// Equality operator: Only considers exact match
    template <typename scene_handle_t, typename intersector_state_t>
    DETRAY_HOST_DEVICE void operator()(state &st,
                                       intersector_state_t &intr_state,
                                       const scene_handle_t & /*geo*/) const {
        if (intr_state.m_is_inside) {
            st.m_pixel.set_color(red);
        } else {
            vector3D dir = vector::normalize(intr_state.m_ray.dir());
            point3D p1{1.0f, 1.0f, 1.0f};
            point3D p2{0.5f, 0.7f, 1.0f};
            const scalar t{0.5f * dir[1] + 1.0f};
            point3D tmp = 255.99f * ((1.0f - t) * p1 + t * p2);

            st.m_pixel.set_color({static_cast<color_depth>(tmp[0]),
                                  static_cast<color_depth>(tmp[1]),
                                  static_cast<color_depth>(tmp[2]),
                                  static_cast<color_depth>(255u)});
        }
    }
};

}  // namespace detray
