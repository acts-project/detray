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
#include "detray/propagator/base_actor.hpp"
#include "detray/tracer/definitions/colors.hpp"
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/detail/material_color_helper.hpp"
#include "detray/tracer/texture/pixel.hpp"

namespace detray {

#if(IS_SOA)
using algebra::storage::operator*;
using algebra::storage::operator/;
using algebra::storage::operator-;
using algebra::storage::operator+;
#endif

/// Calculates the color of a pixel. Starting point of the shader pipeline
struct material_shader : public detray::actor {

    /// Equality operator: Only considers exact match
    template <typename scene_handle_t, typename intersector_state_t>
    DETRAY_HOST_DEVICE void operator()(state &, intersector_state_t &intr_state,
                                       scene_handle_t &sc) const {
        using color_depth = typename decltype(sc.m_pixel)::color_depth;

        if (intr_state.m_is_inside) {
            // auto c = texture::detail::material_color_helper<color_depth>(
            //         intr_state.material());
            auto normal = sc.geometry().mask().normal(
                intr_state.m_intersections[0].local);
            normal = normal + vector3D{1.f, 1.f, 1.f};
            normal = 255.99f * 0.5f * normal;
            sc.m_pixel.set_color({static_cast<color_depth>(normal[0][0]),
                                  static_cast<color_depth>(normal[1][0]),
                                  static_cast<color_depth>(normal[2][0])});
        }
    }
};

}  // namespace detray
