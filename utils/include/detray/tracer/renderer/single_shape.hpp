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
#include "detray/intersection/intersection.hpp"
#include "detray/propagator/base_actor.hpp"

// System include(s)
#include <limits>

namespace detray {

/// Calculates the color of a pixel. Starting point of the shader pipeline
template <typename mask_t>
struct single_shape : detray::actor {

    // dummy output type: will be gone with future update
    using output_type = bool;

    struct state {

        DETRAY_HOST_DEVICE
        state(const mask_t &mask, const point3D &ori, const vector3D &dir,
              scalar min = 0.f,
              scalar max = std::numeric_limits<scalar>::infinity())
            : m_mask{mask}, m_ray{ori, 0.f, dir, 0.f}, m_interval{min, max} {}

        /// The shape to be rendered
        mask_t m_mask;
        /// The input ray into the scene
        detail::ray<transform3D> m_ray;
        /// Interval in which ray intersections are being considered
        std::array<scalar, 2> m_interval;
        /// Resulting intersection
        std::array<line_plane_intersection, 2> m_intersections;
        /// Flag to the obseving colorizer
        bool m_is_inside = false;
    };

    /// Intersect the ray with the mask. The closest intersection is in front of
    /// the @c m_intersections container
    template <typename scene_handle_t>
    DETRAY_HOST_DEVICE void operator()(state &st,
                                       const scene_handle_t & /*geo*/) const {
        vector3D x{1.0f, 0.0f, 0.0f};
        vector3D z{0.0f, 0.0f, 1.0f};
        vector3D t{1.0f, 1.0f, 10.0f};
        transform3D trf{t, z, x};
        st.m_is_inside = place_in_collection(
            st.m_mask.intersector()(st.m_ray, st.m_mask, trf),
            st.m_intersections);
    }

    private:
    template <typename is_container_t>
    DETRAY_HOST_DEVICE bool place_in_collection(
        std::array<typename is_container_t::value_type, 1> &&sfi,
        is_container_t &intersections) const {
        if (sfi[0].status == intersection::status::e_inside) {
            intersections[0] = sfi[0];
            return true;
        } else {
            return false;
        }
    }

    template <typename is_container_t>
    DETRAY_HOST_DEVICE bool place_in_collection(
        std::array<typename is_container_t::value_type, 2> &&solutions,
        is_container_t &intersections) const {
        bool is_valid = false;
        for (auto &sfi : solutions) {
            if (sfi.status == intersection::status::e_inside) {
                intersections.push_back(sfi);
                is_valid = true;
            }
        }
        return is_valid;
    }
};

}  // namespace detray
