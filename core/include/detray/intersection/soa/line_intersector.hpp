/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/line2D.hpp"
#include "detray/definitions/boolean.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"

// System include(s)
#include <type_traits>

namespace detray::soa {

/// A functor to find intersections between trajectory and line mask
template <typename intersection_t>
struct line_intersector {

    /// Linear algebra types
    /// @{
    using algebra = typename intersection_t::algebra;
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

    /// Operator function to find intersections between ray and line mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    //
    /// @return the intersection
    template <typename mask_t, typename surface_t,
              std::enable_if_t<std::is_same_v<typename mask_t::local_frame_type,
                                              line2D<algebra>>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline intersection_t operator()(
        const ray_type &ray, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        intersection_type is;

        // line direction
        const vector3 sz = getter::vector<3>(trf.matrix(), 0u, 2u);

        // line center
        const point3 st = trf.translation();

        // Broadcast ray data
        const auto &pos = ray.pos();
        const auto &dir = ray.dir();
        const vector3 ro{pos[0], pos[1], pos[2]};
        const vector3 rd{dir[0], dir[1], dir[2]};

        // Projection of line to track direction
        const scalar_type zd = vector::dot(sz, rd);

        const scalar_type denom = 1.f - (zd * zd);

        // Case for wire is parallel to track
        if (detray::detail::all_of(denom < 1e-5f)) {
            is.status = decltype(is.status)(false);
            return is;
        }

        // vector from track position to line center
        const vector3 t2l = st - ro;

        // t2l projection on line direction
        const scalar_type t2l_on_line = vector::dot(t2l, sz);

        // t2l projection on track direction
        const scalar_type t2l_on_track = vector::dot(t2l, rd);

        // path length to the point of closest approach on the track
        is.path = (t2l_on_track - t2l_on_line * zd) / denom;

        // point of closest approach on the track
        const point3 m = ro + rd * is.path;
        is.local = mask.to_local_frame(trf, m, rd);
        is.status = mask.is_inside(is.local, mask_tolerance);

        // Early return, in case all intersections are invalid
        if (detray::detail::none_of(is.status)) {
            return is;
        }

        is.sf_desc = sf;

        is.direction = math_ns::signbit(is.path);
        is.volume_link = mask.volume_link();

        // Get incidence angle
        is.cos_incidence_angle = math_ns::abs(zd);

        // Mask the values where the overstepping tolerance was not met
        is.status &= (is.path >= ray.overstep_tolerance());

        return is;
    }

    /// Operator function to find intersections between a ray and a line.
    ///
    /// @tparam mask_t is the input mask type
    ///
    /// @param ray is the input ray trajectory
    /// @param sfi the intersection to be updated
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    template <typename mask_t,
              std::enable_if_t<std::is_same_v<typename mask_t::local_frame_type,
                                              line2D<algebra>>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray, intersection_t &sfi, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {
        sfi = this->operator()(ray, sfi.sf_desc, mask, trf, mask_tolerance);
    }
};

}  // namespace detray::soa
