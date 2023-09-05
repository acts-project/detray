/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"

// System include(s)
#include <type_traits>

namespace detray::soa {

/// A functor to find intersections between straight line and planar surface
template <typename intersection_t>
struct plane_intersector {

    /// linear algebra types
    /// @{
    using transform3_type = typename intersection_t::transform3D;
    using value_type = typename intersection_t::value_t;
    using scalar_type = typename intersection_t::scalar_t;
    using point3 = typename intersection_t::point3D;
    using point2 = typename intersection_t::point2D;
    using vector3 = typename intersection_t::vector3D;
    /// @}

    using intersection_type = intersection_t;
    using ray_type =
        detray::detail::ray<dtransform3D<ALGEBRA_PLUGIN<value_type>>>;

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
    template <typename mask_t, typename surface_t>
    DETRAY_HOST_DEVICE inline intersection_t operator()(
        const ray_type &ray, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        intersection_t is;

        // Retrieve the surface normal & translation (context resolved)
        const auto &sm = trf.matrix();
        // const vector3 sn = getter::vector<3>(sm, 0u, 2u);
        const vector3 sn = sm.z;
        const vector3 st = trf.translation();

        // Broadcast ray data
        const auto &pos = ray.pos();
        const auto &dir = ray.dir();
        const vector3 ro{pos[0], pos[1], pos[2]};
        const vector3 rd{dir[0], dir[1], dir[2]};

        const scalar_type denom = vector::dot(rd, sn);
        const vector3 diff = st - ro;
        is.path = vector::dot(sn, diff) / denom;

        // Check if we divided by zero
        const auto check_sum = is.path.sum();
        if (!std::isnan(check_sum) and !std::isinf(check_sum)) {

            const point3 p3 = ro + is.path * rd;
            is.local = mask.to_local_frame(trf, p3, rd);
            is.status = mask.is_inside(is.local, mask_tolerance);

            // Early return, if no intersection was found
            if (is.status.isEmpty()) {
                return is;
            }

            is.sf_desc = sf;

            is.direction = math_ns::signbit(is.path);
            is.volume_link = mask.volume_link();

            // Get incidene angle
            is.cos_incidence_angle = math_ns::abs(denom);

            // Mask the values where the overstepping tolerance was not met
            is.status &= (is.path >= ray.overstep_tolerance());
        } else {
            is.status = decltype(is.status)(false);
        }

        return is;
    }

    /// Operator function to updtae an intersections between a ray and a planar
    /// surface.
    ///
    /// @tparam mask_t is the input mask type
    ///
    /// @param ray is the input ray trajectory
    /// @param sfi the intersection to be updated
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray, intersection_t &sfi, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {
        sfi = this->operator()(ray, sfi.sf_desc, mask, trf, mask_tolerance);
    }
};

}  // namespace detray::soa
