/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/math.hpp"
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
template <typename T, template <typename> class algebra_t>
struct plain_background : public image_background {

    using transform3D = dtransform3D<algebra_t<T>>;

    template <typename color_depth = unsigned int>
    constexpr texture::color<color_depth> get(
        const detray::ray<transform3D> &) {
        return m_color<color_depth>;
    }

    template <typename color_depth = unsigned int>
    static constexpr auto m_color = texture::white<color_depth>;
};

/// @brief Gradient background as described in
template <typename T, template <typename> class algebra_t>
struct gradient_background : public image_background {

    using point3D = dpoint3D<algebra_t<T>>;
    using vector3D = dvector3D<algebra_t<T>>;
    using transform3D = dtransform3D<algebra_t<T>>;

    template <typename color_depth = unsigned int>
    constexpr texture::color<color_depth> get(
        const detray::ray<transform3D> &ray) {
        vector3D dir = vector::normalize(ray.dir());
        point3D p1{1.0f, 1.0f, 1.0f};
        point3D p2{0.85f, 0.85f, 1.0f};
        const auto t = 0.5f * dir[1] + 1.0f;
        point3D p4 = ((1.0f - t) * p1 + t * p2);
        point3D p3 = 255.99f * p4;

        return {static_cast<color_depth>(p3[0][0]),
                static_cast<color_depth>(p3[1][0]),
                static_cast<color_depth>(p3[2][0]),
                static_cast<color_depth>(255u)};
    }
};

/// @brief Gradient background as described in
template <class image_background_t>
struct inf_plane : public image_background {

    using transform3D = typename image_background_t::transform3D;

    template <typename color_depth = unsigned int>
    constexpr texture::color<color_depth> get(
        const detray::ray<transform3D> &ray) {

#if (IS_SOA)
        if ((ray.dir()[1] < 0.f).isFull()) {
#else
        if (not std::signbit(ray.dir()[1])) {
#endif
            return texture::green<color_depth>;
        } else {
            return image_background_t{}.template get<color_depth>(ray);
        }
    }
};

/// Calculates the color of a pixel. Starting point of the shader pipeline
template <class image_background_t>
struct background_shader : public detray::actor {

    /// Equality operator: Only considers exact match
    template <typename scene_handle_t, typename intersector_state_t>
    DETRAY_HOST_DEVICE void operator()(state &, intersector_state_t &intr_state,
                                       scene_handle_t &sc) const {
        using color_depth = typename decltype(sc.m_pixel)::color_depth;

        if (not intr_state.m_is_inside) {
            sc.m_pixel.set_color(
                image_background_t{}.template get<color_depth>(sc.ray()));
        }
    }
};

}  // namespace detray
