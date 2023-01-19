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
#include "tests/common/tools/intersectors/helix_intersector.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

namespace detail {

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
    template <typename mask_t, typename surface_t>
    DETRAY_HOST_DEVICE inline intersection2D<surface_t, transform3_t>
    operator()(const helix_type &h, surface_t sf, const mask_t &mask,
               const transform3_t &trf,
               const scalar_type mask_tolerance = 0.f) const {

        using intersection_t = intersection2D<surface_t, transform3_t>;
        intersection_t sfi;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{100};
        // Tolerance for convergence
        constexpr scalar_type tol{1e-4f};

        // Get the surface info
        const auto &sm = trf.matrix();
        // Surface normal
        const vector3 sn = getter::vector<3>(sm, 0, 2);
        // Surface translation
        const point3 st = getter::vector<3>(sm, 0, 3);

        // Starting point on the helix for the Newton iteration
        scalar_type s{0.95f * getter::norm(st)};
        scalar_type s_prev{0.9f * s};

        // f(s) = sn * (h.pos(s) - st) == 0
        // Run the iteration on s
        std::size_t n_tries{0};
        while (std::abs(s - s_prev) > tol and n_tries < max_n_tries) {
            // f'(s) = sn * h.dir(s)
            const scalar_type denom{vector::dot(sn, h.dir(s))};
            // No intersection can be found if dividing by zero
            if (denom == 0.) {
                return sfi;
            }
            // x_n+1 = x_n - f(s) / f'(s)
            s_prev = s;
            s -= vector::dot(sn, h.pos(s) - st) / denom;
            ++n_tries;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return sfi;
        }

        // Build intersection struct from helix parameters
        const point3 helix_pos = h.pos(s);

        sfi.path = getter::norm(helix_pos);
        if (sfi.path < h.overstep_tolerance()) {
            return sfi;
        }

        sfi.p3 = helix_pos;
        sfi.p2 = mask.to_local_frame(trf, sfi.p3, h.dir(s));
        sfi.status = mask.is_inside(sfi.p2, mask_tolerance);

        // Compute some additional information if the intersection is valid
        if (sfi.status == intersection::status::e_inside) {
            sfi.surface = sf;
            sfi.direction = std::signbit(vector::dot(st, h.dir(s)))
                                ? intersection::direction::e_opposite
                                : intersection::direction::e_along;
            sfi.volume_link = mask.volume_link();
        }

        return sfi;
    }
};

}  // namespace detail

/// Specialization of the @c helix_intersector for planar surfaces
template <typename transform3_t, typename mask_t>
struct helix_intersector<
    transform3_t, mask_t,
    std::enable_if_t<std::is_same_v<typename mask_t::shape::
                                        template intersector_type<transform3_t>,
                                    plane_intersector<transform3_t>>,
                     void>>
    : public detail::helix_plane_intersector<transform3_t> {

    using intersector_impl = detail::helix_plane_intersector<transform3_t>;

    using scalar_type = typename intersector_impl::scalar_type;
    using matrix_operator = typename intersector_impl::matrix_operator;
    using point3 = typename intersector_impl::point3;
    using vector3 = typename intersector_impl::vector3;
    using helix_type = typename intersector_impl::helix_type;
};

}  // namespace detray