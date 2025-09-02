/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/line2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool resolve_pos>
struct ray_intersector_impl;

/// A functor to find intersections between trajectory and line mask
template <algebra::concepts::aos algebra_t, bool resolve_pos>
struct ray_intersector_impl<line2D<algebra_t>, algebra_t, resolve_pos> {

    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    template <typename surface_descr_t>
    using intersection_type =
        intersection2D<surface_descr_t, algebra_t, resolve_pos>;

    template <typename other_algebra_t>
    using trajectory_type = detail::ray<other_algebra_t>;

    // Maximum number of solutions this intersector can produce
    static constexpr std::uint8_t n_solutions{1u};

    using result_type =
        intersection_point<algebra_t, point3_type, intersection::contains_pos>;

    /// Operator function to find intersections between ray and line mask
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tol negative cutoff for the path
    //
    /// @return the intersection
    DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
        const trajectory_type<algebra_t> &ray, const transform3_type &trf,
        const scalar_type /*overstep_tol*/ = 0.f) const {

        const vector3_type &rd = ray.dir();
        const point3_type &ro = ray.pos();

        // line direction
        const vector3_type &_z = trf.z();
        // line center
        const point3_type &_t = trf.translation();

        // Projection of line to track direction
        const scalar_type zd{vector::dot(_z, rd)};

        const scalar_type denom{1.f - (zd * zd)};

        // Case for wire is parallel to track
        if (denom < 1e-5f) {
            return {};
        }

        // vector from track position to line center
        const point3_type t2l = _t - ro;

        // t2l projection on line direction
        const scalar_type t2l_on_line{vector::dot(t2l, _z)};

        // t2l projection on track direction
        const scalar_type t2l_on_track{vector::dot(t2l, rd)};

        // path length to the point of closest approach on the track
        const scalar_type s{1.f / denom * (t2l_on_track - t2l_on_line * zd)};
        const point3_type glob_pos = ro + s * rd;

        return {s, glob_pos};
    }
};

}  // namespace detray
