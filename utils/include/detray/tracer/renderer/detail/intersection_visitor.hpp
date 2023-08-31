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

namespace detray::detail {

/// Calculates the color of a pixel. Starting point of the shader pipeline
struct intersector_visitor {

    // dummy output type: will be gone with future update
    using output_type = bool;

    template <typename mask_group_t, typename mask_range_t,
              typename is_container_t, typename traj_t, typename surface_t,
              typename transform_container_t>
    DETRAY_HOST_DEVICE inline void operator()(
        const mask_group_t &mask_group, const mask_range_t &mask_range,
        std::array<line_plane_intersection, 2> &inersections,
        const detail::ray<transform3D> &ray, const surface_t &surface,
        const transform3D &ctf, const std::array<scalar, 2> /*interval*/,
        const scalar mask_tolerance = 0.f) const {
        // Run over the masks that belong to the surface (only one can be hit)
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {
            if (place_in_collection(
                    mask.intersector()(traj, mask, ctf, mask_tolerance),
                    is_container)) {
                return;
            }
        }
    }

    private:
    template <typename is_container_t>
    DETRAY_HOST_DEVICE bool place_in_collection(
        typename is_container_t::value_type &&sfi,
        is_container_t &intersections) const {
        if (sfi.status) {
            intersections.push_back(sfi);
        }
        return sfi.status;
    }

    template <typename is_container_t>
    DETRAY_HOST_DEVICE bool place_in_collection(
        std::array<typename is_container_t::value_type, 2> &&solutions,
        is_container_t &intersections) const {
        bool is_valid = false;
        for (auto &sfi : solutions) {
            if (sfi.status) {
                intersections.push_back(sfi);
                is_valid = true;
            }
        }
        return is_valid;
    }
};

}  // namespace detray::detail
