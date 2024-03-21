/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/coordinates/line2D.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/intersection/intersection.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, typename algebra_t>
struct ray_intersector_impl;

/// A functor to find intersections between trajectory and line mask
template <typename algebra_t>
struct ray_intersector_impl<line2D<algebra_t>, algebra_t> {

    using transform3_type = algebra_t;
    using scalar_type = typename transform3_type::scalar_type;
    using point3 = typename transform3_type::point3;
    using vector3 = typename transform3_type::vector3;

    template <typename surface_descr_t>
    using intersection_type = intersection2D<surface_descr_t, algebra_t>;
    using ray_type = detail::ray<algebra_t>;

    /// Operator function to find intersections between ray and line mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_descr_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tol negative cutoff for the path
    //
    /// @return the intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const ray_type &ray, const surface_descr_t &sf, const mask_t &mask,
        const transform3_type &trf, const scalar_type mask_tolerance = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;

        // line direction
        const vector3 _z = getter::vector<3>(trf.matrix(), 0u, 2u);

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
        if (is.path >= overstep_tol) {

            // point of closest approach on the track
            const point3 m = _p + _d * A;

            is.local = mask.to_local_frame(trf, m, _d);

            is.status = mask.is_inside(is.local, mask_tolerance);

            // prepare some additional information in case the intersection
            // is valid
            if (is.status == intersection::status::e_inside) {
                is.sf_desc = sf;

                is.direction = detail::signbit(is.path)
                                   ? intersection::direction::e_opposite
                                   : intersection::direction::e_along;
                is.volume_link = mask.volume_link();

                // Get incidence angle
                is.cos_incidence_angle = math::abs(zd);
            }
        }
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
    /// @param overstep_tol negative cutoff for the path
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray, intersection_type<surface_descr_t> &sfi,
        const mask_t &mask, const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        sfi = this->operator()(ray, sfi.sf_desc, mask, trf, mask_tolerance,
                               overstep_tol);
    }
};

}  // namespace detray
