/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical2.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

/// @brief A functor to find intersections between a straight line and a
/// cylindrical portal surface.
///
/// With the way the navigation works, only the closest one of the two possible
/// intersection points is needed in the case of a cylinderical portal surface.
template <typename transform3_t>
struct cylinder_portal_intersector : public cylinder_intersector<transform3_t> {

    /// linear algebra types
    /// @{
    using scalar_type = typename transform3_t::scalar_type;
    using point3 = typename transform3_t::point3;
    using point2 = typename transform3_t::point2;
    using vector3 = typename transform3_t::vector3;
    /// @}
    using ray_type = detail::ray<transform3_t>;

    /// Operator function to find intersections between ray and cylinder mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam transform_t is the input transform type
    ///
    /// @param ray is the input ray trajectory
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tolerance is the tolerance for track overstepping
    ///
    /// @return the closest intersection
    template <
        typename mask_t, typename surface_t,
        std::enable_if_t<std::is_same_v<typename mask_t::measurement_frame_type,
                                        cylindrical2<transform3_t>>,
                         bool> = true>
    DETRAY_HOST_DEVICE inline std::array<
        line_plane_intersection<surface_t, transform3_t>, 1>
    operator()(const ray_type &ray, const surface_t sf, const mask_t &mask,
               const transform3_t &trf,
               const scalar_type mask_tolerance = 0.f) const {

        using intersection_t = line_plane_intersection<surface_t, transform3_t>;
        std::array<intersection_t, 1> ret;

        // Intersecting the cylinder from the inside yield one intersection
        // along the direction of the track and one behind it
        const auto qe = this->solve_intersection(ray, mask, trf);

        // Find the clostest valid intersection
        if (qe.solutions() > 0) {
            // Only the closest intersection that is outside the overstepping
            // tolerance is needed
            const scalar_type t{(qe.smaller() > ray.overstep_tolerance())
                                    ? qe.smaller()
                                    : qe.larger()};
            ret[0] = this->template build_candidate<intersection_t>(
                ray, mask, trf, t, mask_tolerance);
            ret[0].surface = sf;
        } else {
            ret[0].status = intersection::status::e_missed;
        }

        return ret;
    }
};

}  // namespace detray