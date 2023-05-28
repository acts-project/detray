/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/propagator/actor_chain.hpp"
#include "detray/tracer/renderer/intersector.hpp"
#include "detray/tracer/renderer/raw_image.hpp"
#include "detray/tracer/renderer/single_shape.hpp"
#include "detray/tracer/texture/pixel.hpp"

namespace detray {

/// Executes the rendering steps sequentially
template <typename intersector_t, typename background_shader_t, typename... shaders_t>
using rendering_pipeline = actor_chain<dtuple, composite_actor<dtuple, intersector_t, background_shader_t, shaders_t...>>;

/// State that is passed through the pipeline per ray
///
/// Contains a pointer to the geometry and the ray that this pipeline instance
/// renders.
struct scene_handle {

    struct config {};

    template <typename geometry_t, typename color_depth, typename pixel_coord>
    struct state {

        DETRAY_HOST_DEVICE
        state(const geometry_t &geo, 
              const raw_image<color_depth> &im, 
              const detail::ray<transform3D> &ray, 
              const pixel_coord x, const pixel_coord y)
            : m_geo{&geo}, m_image{&im}, m_ray{&ray}, m_pixel{{x, y}} {}

        /// Threadsafe interface
        /// @{
        const detail::ray<transform3D> &ray() const { return *m_ray; }
        const geometry_t &geometry() const { return *m_geo; }
        /// @}

        /// The geometry handle
        const geometry_t *m_geo;
        /// The image handle
        const raw_image<color_depth> *m_image;
        /// The ray handle
        const detail::ray<transform3D> *m_ray;
        /// The pixel for this ray
        texture::pixel<pixel_coord, color_depth> m_pixel;
    };

#if __clang__
    template <typename geometry_t>
    DETRAY_HOST_DEVICE state(const geometry_t &geo,
                             const detail::ray<transform3D> &ray)
        ->state<geometry_t>;
#endif
};

}  // namespace detray
