/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <cmath>
#include <type_traits>

#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/utils/quadratic_equation.hpp"

namespace detray {

namespace detail {

struct unbound;

}

/// @brief Intersection implementation for cylinder surfaces using ray
/// trajectories.
template <template <typename, std::size_t> class array_t = darray>
struct ray_concentric_cylinder_intersector {

    using intersection_type = line_plane_intersection;

    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;
    using point2 = __plugin::point2<detray::scalar>;
    using cylindrical2 = __plugin::cylindrical2<detray::scalar>;

    /// Intersection method for cylindrical surfaces
    ///
    /// @tparam transform_t The type of placement matrix of the cylinder surface
    /// @tparam track_t The type of the track caryying also the context object
    /// @tparam mask_t The mask type applied to the local frame
    ///
    /// Contextual part:
    /// @param trf the transform of the surface surface to be intersected @note
    /// is ignored
    /// @param track the track information
    ///
    /// Non-contextual part:
    /// @param mask the local mask
    /// @param tolerance is the mask specific tolerance
    ///
    /// @return the intersection with optional parameters
    template <
        typename transform_t, typename track_t, typename mask_t,
        std::enable_if_t<
            std::is_same_v<typename mask_t::local_type, cylindrical2> or
                std::is_same_v<typename mask_t::local_type, detail::unbound>,
            bool> = true,
        std::enable_if_t<not std::is_same_v<track_t, detail::ray>, bool> = true>
    DETRAY_HOST_DEVICE inline intersection_type intersect(
        const transform_t &trf, const track_t &track, const mask_t &mask,
        const typename mask_t::mask_tolerance &tolerance =
            mask_t::within_epsilon) const {

        return intersect(trf, detail::ray(track), mask, tolerance,
                         track.overstep_tolerance());
    }

    /// Intersection method for cylindrical surfaces
    ///
    /// @tparam transform_t The type of placement matrix of the cylinder surface
    /// @tparam mask_t The mask type applied to the local frame
    ///
    /// Contextual part:
    /// @param trf the transform of the surface to be intersected
    /// @param ro the origin of the ray
    /// @param rd the direction of the ray
    ///
    /// Non-contextual part:
    /// @param mask the local mask
    /// @param tolerance is the mask specific tolerance
    /// @param overstep_tolerance  is the stepping specific tolerance
    ///
    /// @return the intersection with optional parameters
    template <
        typename transform_t, typename mask_t,
        std::enable_if_t<
            std::is_same_v<typename mask_t::local_type, cylindrical2> or
                std::is_same_v<typename mask_t::local_type, detail::unbound>,
            bool> = true>
    DETRAY_HOST_DEVICE inline intersection_type intersect(
        const transform_t & /*trf*/, const detail::ray &ray, const mask_t &mask,
        const dindex /*volume_index*/ = dindex_invalid,
        const typename mask_t::mask_tolerance & /*tolerance*/ =
            mask_t::within_epsilon,
        const scalar overstep_tolerance = 0.) const {

        using local_frame = typename mask_t::local_type;

        const scalar r{mask[0]};

        // Two points on the line, thes are in the cylinder frame
        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();
        const point3 &l0 = ro;
        const point3 l1 = ro + rd;

        // swap coorinates x/y for numerical stability
        const bool swap_x_y = std::abs(rd[0]) < scalar{1e-3};

        std::size_t _x = swap_x_y ? 1 : 0;
        std::size_t _y = swap_x_y ? 0 : 1;
        const scalar k{(l0[_y] - l1[_y]) / (l0[_x] - l1[_x])};
        const scalar d{l1[_y] - k * l1[_x]};

        quadratic_equation<scalar> qe = {(1 + k * k), 2 * k * d, d * d - r * r};
        auto qe_solution = qe();

        if (std::get<0>(qe_solution) > overstep_tolerance) {
            array_t<point3, 2> candidates;
            const auto u01 = std::get<1>(qe_solution);
            array_t<scalar, 2> t01 = {0., 0.};

            candidates[0][_x] = u01[0];
            candidates[0][_y] = k * u01[0] + d;
            t01[0] = (candidates[0][_x] - ro[_x]) / rd[_x];
            candidates[0][2] = ro[2] + t01[0] * rd[2];

            candidates[1][_x] = u01[1];
            candidates[1][_y] = k * u01[1] + d;
            t01[1] = (candidates[1][_x] - ro[_x]) / rd[_x];
            candidates[1][2] = ro[2] + t01[1] * rd[2];

            // Chose the index, take the smaller positive one
            const std::size_t cindex =
                (t01[0] < t01[1] and t01[0] > overstep_tolerance)
                    ? 0
                    : (t01[0] < overstep_tolerance and
                               t01[1] > overstep_tolerance
                           ? 1
                           : 0);
            if (t01[0] > overstep_tolerance or t01[1] > overstep_tolerance) {
                intersection_type is;
                is.p3 = candidates[cindex];
                is.path = t01[cindex];

                is.p2 = point2{r * getter::phi(is.p3), is.p3[2]};
                is.status = mask.template is_inside<local_frame>(is.p3);
                const scalar rdir{getter::perp(is.p3 + scalar{0.1} * rd)};
                is.direction = rdir > r ? intersection::direction::e_along
                                        : intersection::direction::e_opposite;
                is.link = mask.volume_link();
                return is;
            }
        }
        return intersection_type{};
    }
};

}  // namespace detray
