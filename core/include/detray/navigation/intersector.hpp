/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/navigation/detail/trajectories.hpp"
#include "detray/navigation/intersection/helix_intersector.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"

namespace detray {

/// @brief Intersection interface for detector surfaces.
///
/// Composes the different intersector options into a unifyed interface
template <typename algebra_t, typename shape_t>
struct intersector {

    using transform3_type = algebra_t;
    using scalar_type = typename transform3_type::scalar_type;

    /// How to intersect surfaces with rays
    using ray_intersector_type = ray_intersector<algebra_t, shape_t>;

    /// How to intersect surfaces with helices
    using helix_intersector_type = helix_intersector<algebra_t, shape_t>;

    /// @returns the intersection(s) between a surface and the ray @param ray
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline decltype(auto) operator()(
        const detail::ray<algebra_t> &ray, const surface_descr_t &sf,
        const mask_t &mask, const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        return ray_intersector_type{}(ray, sf, mask, trf, mask_tolerance,
                                      overstep_tol);
    }

    /// @returns the intersection(s) between a surface and the helix @param h
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline decltype(auto) operator()(
        const detail::helix<algebra_t> &h, const surface_descr_t &sf,
        const mask_t &mask, const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f, const scalar_type = 0.f) const {

        return helix_intersector_type{}(h, sf, mask, trf, mask_tolerance);
    }
};

}  // namespace detray
