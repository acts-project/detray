/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/propagator/actor_chain.hpp"
#include "detray/tracer/renderer/colorizer.hpp"
#include "detray/tracer/renderer/intersector.hpp"
#include "detray/tracer/renderer/single_shape.hpp"

namespace detray {

/// Executes the rendering steps sequentially
template <typename... shaders_t>
using rendering_pipeline = actor_chain<dtuple, shaders_t...>;

/// State that is passed through the pipeline per ray
///
/// Contains a pointer to the geometry and the ray that this pipeline instance
/// renders.
struct scene_handle {

    struct config {};

    template <typename geometry_t>
    struct state {

        DETRAY_HOST_DEVICE
        state(const geometry_t &geo, const detail::ray<transform3D> &ray)
            : m_geo{&geo}, m_ray{&ray} {}

        /// Threadsafe interface
        /// @{
        const detail::ray<transform3D> &ray() const { return *m_ray; }
        const geometry_t &geometry() const { return *m_geo; }
        /// @}

        /// The geometry handle
        const geometry_t *m_geo;
        /// The ray handle
        const detail::ray<transform3D> *m_ray;
    };

#if __clang__
    template <typename geometry_t>
    DETRAY_HOST_DEVICE state(const geometry_t &geo,
                             const detail::ray<transform3D> &ray)
        ->state<geometry_t>;
#endif
};

}  // namespace detray
