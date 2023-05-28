/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/tracer/definitions/colors.hpp"
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/pixel.hpp"

namespace detray {

/// @brief Image background class tag
struct image_background {};

/// @brief Single color image background
template <typename color_depth = uint8_t>
struct plain_background : public image_background {

    using color_t = texture::color<color_depth>;

    constexpr color_t operator()(const detray::detail::ray<transform3D> &) {
        return m_color;
    }

    color_t m_color = texture::white<color_depth>;
};

/// @brief Gradient background as described in
template <typename color_depth = uint8_t>
struct gradient_background : public image_background {

    using color_t = texture::color<color_depth>;

    constexpr color_t operator()(const detray::detail::ray<transform3D> &ray) {
        vector3D dir = vector::normalize(ray.dir());
        point3D p1{1.0f, 1.0f, 1.0f};
        point3D p2{0.85f, 0.85f, 1.0f};
        const scalar t{0.5f * dir[1] + 1.0f};
        point3D p3 = 255.99f * ((1.0f - t) * p1 + t * p2);

        return {
            static_cast<color_depth>(p3[0]), static_cast<color_depth>(p3[1]),
            static_cast<color_depth>(p3[2]), static_cast<color_depth>(255u)};
    }
};

/// Calculates the color of a pixel. Starting point of the shader pipeline
template <template <typename> class image_background_t = plain_background>
struct background_shader : public detray::actor {

    /// Equality operator: Only considers exact match
    template <typename scene_handle_t, typename intersector_state_t>
    DETRAY_HOST_DEVICE void operator()(state &,
                                       intersector_state_t &intr_state,
                                       scene_handle_t &sc) const {
        using color_depth = typename decltype(sc.m_pixel)::color_depth;

        if (not intr_state.m_is_inside) {
            sc.m_pixel.set_color(image_background_t<color_depth>{}(sc.ray()));
        }
    }
};

}  // namespace detray
