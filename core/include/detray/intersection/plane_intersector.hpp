/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

/// A functor to find intersections between straight line and planar surface
template <typename transform3_t>
struct plane_intersector {

    /// linear algebra types
    /// @{
    using scalar_type = typename transform3_t::scalar_type;
    using point3 = typename transform3_t::point3;
    using point2 = typename transform3_t::point2;
    using vector3 = typename transform3_t::vector3;
    /// @}
    using ray_type = detail::ray<transform3_t>;

    /// Operator function to find intersections between ray and planar mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    ///
    /// @return the intersection
    template <
        typename mask_t, typename surface_t,
        std::enable_if_t<std::is_same_v<typename mask_t::loc_point_t, point2>,
                         bool> = true>
    DETRAY_HOST_DEVICE inline intersection2D<surface_t, transform3_t>
    operator()(const ray_type &ray, const surface_t sf, const mask_t &mask,
               const transform3_t &trf,
               const scalar_type mask_tolerance = 0.f) const {

        using intersection_t = intersection2D<surface_t, transform3_t>;
        intersection_t is;

        // Retrieve the surface normal & translation (context resolved)
        const auto &sm = trf.matrix();
        const vector3 sn = getter::vector<3>(sm, 0, 2);
        const vector3 st = getter::vector<3>(sm, 0, 3);

        // Intersection code
        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();
        const scalar_type denom = vector::dot(rd, sn);
        if (denom != 0.f) {
            is.path = vector::dot(sn, st - ro) / denom;

            // Intersection is not valid for navigation - return early
            if (is.path < ray.overstep_tolerance()) {
                return is;
            }

            is.surface = sf;
            is.p3 = ro + is.path * rd;
            is.p2 = mask.to_local_frame(trf, is.p3, ray.dir());
            is.status = mask.is_inside(is.p2, mask_tolerance);

            // prepare some additional information in case the intersection
            // is valid
            if (is.status == intersection::status::e_inside) {
                is.direction = std::signbit(is.path)
                                   ? intersection::direction::e_opposite
                                   : intersection::direction::e_along;
                is.volume_link = mask.volume_link();

                // Get incidene angle
                is.cos_incidence_angle = std::abs(denom);
            }
        } else {
            is.status = intersection::status::e_missed;
        }

        return is;
    }

    /// Operator function to updtae an intersections between a ray and a planar
    /// surface.
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sfi the intersection to be updated
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    template <
        typename mask_t, typename surface_t,
        std::enable_if_t<std::is_same_v<typename mask_t::loc_point_t, point2>,
                         bool> = true>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray, intersection2D<surface_t, transform3_t> &sfi,
        const mask_t &mask, const transform3_t &trf,
        const scalar_type mask_tolerance = 0.f) const {
        sfi = this->operator()(ray, sfi.surface, mask, trf, mask_tolerance);
    }
};

}  // namespace detray