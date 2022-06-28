/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
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

/** A functor to find intersections between trajectory and planar mask
 */
struct plane_intersector {

    using intersection_type = line_plane_intersection;
    using output_type = std::array<intersection_type, 1>;
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /** Operator function to find intersections between ray and planar mask
     *
     * @tparam mask_t is the input mask type
     * @tparam transform_t is the input transform type
     *
     * @param ray is the input ray trajectory
     * @param mask is the input mask
     * @param trf is the transform
     * @param edge_tolerance is the tolerance for mask edges
     * @param overstep_tolerance is the tolerance for track overstepping
     *
     * @return the intersection
     */
    template <typename mask_t, typename transform_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const detail::ray &ray, const mask_t &mask, const transform_t &trf,
        const scalar edge_tolerance = 0,
        const scalar overstep_tolerance = 0.) const {

        output_type ret;

        using local_frame = typename mask_t::local_type;

        // Retrieve the surface normal & translation (context resolved)
        const auto &sm = trf.matrix();
        const vector3 sn = getter::vector<3>(sm, 0, 2);
        const vector3 st = getter::vector<3>(sm, 0, 3);

        // Intersection code
        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();
        const scalar denom = vector::dot(rd, sn);
        if (denom != scalar{0.}) {
            intersection_type &is = ret[0];
            is.path = vector::dot(sn, st - ro) / denom;
            is.p3 = ro + is.path * rd;
            constexpr local_frame local_converter{};
            is.p2 = local_converter(trf, is.p3);
            is.status =
                mask.template is_inside<local_frame>(is.p2, edge_tolerance);
            is.direction = is.path > overstep_tolerance
                               ? intersection::direction::e_along
                               : intersection::direction::e_opposite;
            is.link = mask.volume_link();
        }
        return ret;
    }

    /** Operator function to find intersections between helix and planar mask
     *
     * @tparam mask_t is the input mask type
     * @tparam transform_t is the input transform type
     *
     * @param h is the input helix trajectory
     * @param mask is the input mask
     * @param trf is the transform
     * @param edge_tolerance is the tolerance for mask edges
     * @param overstep_tolerance is the tolerance for track overstepping
     *
     * @return the intersection
     */
    template <typename mask_t, typename transform_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const detail::helix &h, const mask_t &mask, const transform_t &trf,
        const scalar edge_tolerance = 0,
        const scalar overstep_tolerance = 0.) const {

        (void)overstep_tolerance;

        using local_frame = typename mask_t::local_type;

        output_type ret;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{100};
        // Tolerance for convergence
        constexpr scalar tol{1e-3};

        // Get the surface info
        const auto &sm = trf.matrix();
        // Surface normal
        const vector3 sn = getter::vector<3>(sm, 0, 2);
        // Surface translation
        const point3 st = getter::vector<3>(sm, 0, 3);

        // Starting point on the helix for the Newton iteration
        scalar s{getter::norm(sn) - scalar{0.1}};
        scalar s_prev{s - scalar{0.1}};

        // f(s) = sn * (h.pos(s) - st) == 0
        // Run the iteration on s
        std::size_t n_tries{0};
        while (std::abs(s - s_prev) > tol and n_tries < max_n_tries) {
            // f'(s) = sn * h.dir(s)
            const scalar denom{vector::dot(sn, h.dir(s))};
            // No intersection can be found if dividing by zero
            if (denom == 0.) {
                return ret;
            }
            // x_n+1 = x_n - f(s) / f'(s)
            s_prev = s;
            s -= vector::dot(sn, h.pos(s) - st) / denom;
            ++n_tries;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return ret;
        }

        // Build intersection struct from helix parameter s
        intersection_type &is = ret[0];
        const point3 helix_pos = h.pos(s);

        is.path = getter::norm(helix_pos);
        is.p3 = helix_pos;
        constexpr local_frame local_converter{};
        is.p2 = local_converter(trf, is.p3);

        is.status = mask.template is_inside<local_frame>(is.p2, edge_tolerance);
        is.direction = vector::dot(st, h.dir(s)) > scalar{0.}
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();

        return ret;
    }
};

}  // namespace detray