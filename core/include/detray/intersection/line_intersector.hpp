/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/line2.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

/// A functor to find intersections between trajectory and line mask
template <typename transform3_t>
struct line_intersector {

    /// Transformation matching this struct
    using scalar_type = typename transform3_t::scalar_type;
    using point3 = typename transform3_t::point3;
    using vector3 = typename transform3_t::vector3;
    using ray_type = detail::ray<transform3_t>;

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
    template <
        typename mask_t, typename surface_t,
        std::enable_if_t<std::is_same_v<typename mask_t::measurement_frame_type,
                                        line2<transform3_t>>,
                         bool> = true>
    DETRAY_HOST_DEVICE inline line_plane_intersection<surface_t, transform3_t>
    operator()(const ray_type &ray, const surface_t sf, const mask_t &mask,
               const transform3_t &trf,
               const scalar_type mask_tolerance = 0.f) const {

        using intersection_t = line_plane_intersection<surface_t, transform3_t>;
        intersection_t is;

        // line direction
        const vector3 _z = getter::vector<3>(trf.matrix(), 0, 2);

        // line center
        const point3 _t = trf.translation();

        // track direction
        const vector3 _d = ray.dir();

        // track position
        const point3 _p = ray.pos();

        // Projection of line to track direction
        const scalar_type zd{vector::dot(_z, _d)};

        const scalar_type denom{1.f - (zd * zd)};

        // Case for wire is parallel to track
        if (denom < 1e-5f) {
            is.status = intersection::status::e_missed;
            return is;
        }

        // vector from track position to line center
        const auto t2l = _t - _p;

        // t2l projection on line direction
        const scalar_type t2l_on_line{vector::dot(t2l, _z)};

        // t2l projection on track direction
        const scalar_type t2l_on_track{vector::dot(t2l, _d)};

        // path length to the point of closest approach on the track
        const scalar_type A{1.f / denom * (t2l_on_track - t2l_on_line * zd)};

        is.path = A;
        // Intersection is not valid for navigation - return early
        if (is.path < ray.overstep_tolerance()) {
            return is;
        }

        // distance to the point of closest approarch on the
        // line from line center
        const scalar_type B{zd * A - t2l_on_line};

        // point of closest approach on the track
        const vector3 m = _p + _d * A;
        is.p3 = m;

        // For the radial cross section, the calculation does not need to be
        // repeated
        if constexpr (mask_t::shape::square_cross_sect) {
            // This is cartesian3 for square cross section
            auto loc3D = mask.to_local_frame(trf, is.p3);
            is.status = mask.is_inside(loc3D, mask_tolerance);

            // Determine the measurement point from the local point

            // assign the sign depending on the position w.r.t line surface
            // Right: -1
            // Left: 1
            const auto r = vector::cross(_z, _d);
            is.p2[0] = -std::copysign(getter::perp(loc3D), vector::dot(r, t2l));
        } else {
            // local frame and measurement frame are identical
            is.p2 = mask.to_measurement_frame(trf, is.p3, _d);
            is.status = mask.is_inside(is.p2, mask_tolerance);
        }

        // prepare some additional information in case the intersection
        // is valid
        if (is.status == intersection::status::e_inside) {
            is.p2[1] = B;
            is.surface = sf;
            is.direction = std::signbit(is.path)
                               ? intersection::direction::e_opposite
                               : intersection::direction::e_along;
            is.volume_link = mask.volume_link();

            // Get incidence angle
            is.cos_incidence_angle = std::abs(zd);
        }

        return is;
    }

    /// Operator function to find intersections between a ray and a line.
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
        std::enable_if_t<std::is_same_v<typename mask_t::measurement_frame_type,
                                        line2<transform3_t>>,
                         bool> = true>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray,
        line_plane_intersection<surface_t, transform3_t> &sfi,
        const mask_t &mask, const transform3_t &trf,
        const scalar_type mask_tolerance = 0.f) const {
        sfi = this->operator()(ray, sfi.surface, mask, trf, mask_tolerance);
    }
};

}  // namespace detray