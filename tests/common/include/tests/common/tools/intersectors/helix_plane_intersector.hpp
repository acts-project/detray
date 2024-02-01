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
#include "detray/intersection/plane_intersector.hpp"
#include "tests/common/tools/intersectors/helix_intersector.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <type_traits>

namespace detray {

namespace detail {

/// @brief Intersection implementation for helical trajectories with planar
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <typename intersection_t>
struct helix_plane_intersector {

    using transform3_type = typename intersection_t::transform3_type;
    using scalar_type = typename transform3_type::scalar_type;
    using matrix_operator = typename transform3_type::matrix_actor;
    using point3 = typename transform3_type::point3;
    using vector3 = typename transform3_type::vector3;
    using helix_type = detail::helix<transform3_type>;

    /// Operator function to find intersections between helix and planar mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_desc_t is the input surface descriptor type
    ///
    /// @param h is the input helix trajectory
    /// @param sf_desc is the surface descriptor
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tolerance is the tolerance for track overstepping
    ///
    /// @return the intersection
    template <typename mask_t, typename surface_desc_t>
    DETRAY_HOST_DEVICE inline intersection_t operator()(
        const helix_type &h, const surface_desc_t &sf_desc, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        intersection_t sfi;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{1000u};

        // Get the surface info
        const auto &sm = trf.matrix();
        // Surface normal
        const vector3 sn = getter::vector<3>(sm, 0u, 2u);
        // Surface translation
        const point3 st = getter::vector<3>(sm, 0u, 3u);

        const point3 &ro = h._pos;
        const vector3 &rd = h._t0;

        // Starting point on the helix for the Newton iteration
        scalar_type s = vector::dot(sn, st - ro) / vector::dot(rd, sn);
        scalar_type s_min, s_max;
        if (s < 0.f) {
            s_min = 2 * s;
            s_max = 0.f;
        } else {
            s_min = 0.f;
            s_max = 2 * s;
        }
        scalar_type ds_prev = math::abs(s_max - s_min);
        scalar_type ds = ds_prev;
        scalar_type f, dfds;

        // Bisection + Newton-Raphson method
        // Reference: 460-461 page of [Numerical Recipes. 3rd edition]
        for (std::size_t i = 0u; i < max_n_tries + 1; i++) {

            if (i == max_n_tries) {
                return sfi;
            }

            if (((s - s_max) * dfds - f) * ((s - s_min) * dfds - f) > 0.f ||
                math::abs(2.f * ds) > math::abs(ds_prev)) {
                ds_prev = ds;
                ds = 0.5f * (s_max - s_min);
                s = s_min + ds;

                if (math::abs(s - s_min) <= convergence_tolerance) {
                    break;
                }
            } else {
                ds_prev = ds;
                ds = f / dfds;
                scalar_type tmp = s;
                s -= ds;

                if (math::abs(tmp - s) <= convergence_tolerance) {
                    break;
                }
            }

            if (math::abs(ds) < convergence_tolerance) {
                break;  // Convergence criterion.
            }
            // f'(s) = sn * h.dir(s)
            dfds = vector::dot(sn, h.dir(s));
            // No intersection can be found if dividing by zero
            if (dfds == 0.f) {
                return sfi;
            }

            // f(s) = sn * (h.pos(s) - st) == 0
            f = vector::dot(sn, h.pos(s) - st);

            // The one new function evaluation per iteration.
            if (f < 0.f) {
                s_min = s;
            } else {
                s_max = s;
            }
        }

        // Build intersection struct from helix parameters
        sfi.path = s;
        sfi.local = mask.to_local_frame(trf, h.pos(s), h.dir(s));
        sfi.status = mask.is_inside(sfi.local, mask_tolerance);

        // Compute some additional information if the intersection is valid
        if (sfi.status == intersection::status::e_inside) {
            sfi.sf_desc = sf_desc;
            sfi.direction = math::signbit(s)
                                ? intersection::direction::e_opposite
                                : intersection::direction::e_along;
            sfi.volume_link = mask.volume_link();
        }

        return sfi;
    }

    /// Tolerance for convergence
    scalar_type convergence_tolerance{1e-3f};
};

}  // namespace detail

/// Specialization of the @c helix_intersector for planar surfaces
template <typename intersection_t, typename mask_t>
struct helix_intersector<
    intersection_t, mask_t,
    std::enable_if_t<
        std::is_same_v<
            typename mask_t::shape::template intersector_type<intersection_t>,
            plane_intersector<intersection_t>>,
        void>> : public detail::helix_plane_intersector<intersection_t> {

    using intersector_impl = detail::helix_plane_intersector<intersection_t>;

    using scalar_type = typename intersector_impl::scalar_type;
    using matrix_operator = typename intersector_impl::matrix_operator;
    using point3 = typename intersector_impl::point3;
    using vector3 = typename intersector_impl::vector3;
    using helix_type = typename intersector_impl::helix_type;
};

}  // namespace detray
