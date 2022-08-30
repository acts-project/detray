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
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

/// @brief Intersection implementation for helical trajectories with planar
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <typename transform3_t>
struct helix_plane_intersector {

    using scalar_type = typename transform3_t::scalar_type;
    using matrix_operator = typename transform3_t::matrix_actor;
    using point3 = typename transform3_t::point3;
    using vector3 = typename transform3_t::vector3;
    using helix_type = detail::helix<transform3_t>;
    using intersection_type = line_plane_intersection;
    using output_type = std::array<intersection_type, 2>;

    /// Operator function to find intersections between helix and planar mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam transform_t is the input transform type
    ///
    /// @param h is the input helix trajectory
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tolerance is the tolerance for track overstepping
    ///
    /// @return the intersection
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const helix_type &h, const mask_t &mask, const transform3_t &trf,
        const scalar mask_tolerance = 0) const {

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
        is.p2 = local_converter.global_to_local(trf, is.p3, h.dir(s));

        is.status = mask.template is_inside<local_frame>(is.p2, mask_tolerance);
        is.direction = vector::dot(st, h.dir(s)) > scalar{0.}
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();

        return ret;
    }
};

}  // namespace detray