/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical2.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "tests/common/tools/intersectors/helix_intersector.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <type_traits>

namespace detray {

namespace detail {

/// @brief Intersection implementation for cylinder surfaces using helical
/// trajectories.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <typename transform3_t>
struct helix_cylinder_intersector {

    using scalar_type = typename transform3_t::scalar_type;
    using matrix_operator = typename transform3_t::matrix_actor;
    using point2 = typename transform3_t::point2;
    using point3 = typename transform3_t::point3;
    using vector3 = typename transform3_t::vector3;
    using helix_type = detail::helix<transform3_t>;

    /// Operator function to find intersections between helix and cylinder mask
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
    DETRAY_HOST_DEVICE inline std::array<
        line_plane_intersection<surface_t, transform3_t>, 2>
    operator()(const helix_type &h, surface_t sf, const mask_t &mask,
               const transform3_t &trf,
               const scalar_type mask_tolerance = 0) const {

        using intersection_t = line_plane_intersection<surface_t, transform3_t>;
        std::array<intersection_t, 2> ret;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{100};
        // Tolerance for convergence
        constexpr scalar_type tol{1e-3f};

        // Get the surface placement
        const auto &sm = trf.matrix();
        // Cylinder z axis
        const vector3 sz = getter::vector<3>(sm, 0, 2);
        // Cylinder centre
        const point3 sc = getter::vector<3>(sm, 0, 3);

        // Starting point on the helix for the Newton iteration
        // The mask is a cylinder -> it provides its radius as the first value
        const scalar_type r{mask[cylinder2D<>::e_r]};
        // Helix path length parameter
        scalar_type s{r * getter::perp(h.dir(tol))};
        // Path length in the previous iteration step
        scalar_type s_prev{s - 0.1f};

        // f(s) = ((h.pos(s) - sc) x sz)^2 - r^2 == 0
        // Run the iteration on s
        std::size_t n_tries{0};
        while (std::abs(s - s_prev) > tol and n_tries < max_n_tries) {

            // f'(s) = 2 * ( (h.pos(s) - sc) x sz) * (h.dir(s) x sz) )
            const vector3 crp = vector::cross(h.pos(s) - sc, sz);
            const scalar_type denom{
                2.f * vector::dot(crp, vector::cross(h.dir(s), sz))};
            // No intersection can be found if dividing by zero
            if (denom == 0.f) {
                return ret;
            }
            // x_n+1 = x_n - f(s) / f'(s)
            s_prev = s;
            s -= (vector::dot(crp, crp) - r * r) / denom;

            ++n_tries;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return ret;
        }

        // Build intersection struct from helix parameter s
        intersection_t &is = ret[0];
        const point3 helix_pos = h.pos(s);

        is.path = getter::norm(helix_pos);
        if (is.path < h.overstep_tolerance()) {
            return ret;
        }
        is.p3 = helix_pos;
        is.p2 = mask.to_local_frame(trf, is.p3, h.dir(s));
        is.status = mask.is_inside(is.p2, mask_tolerance);

        // Additionally check radial position for Newton solution
        const scalar_type radial_pos{getter::perp(trf.point_to_local(is.p3))};
        const bool r_check =
            std::abs(r - radial_pos) <
            mask_tolerance + 5.f * std::numeric_limits<scalar_type>::epsilon();
        if (not r_check) {
            is.status = intersection::status::e_outside;
        }

        // Compute some additional information if the intersection is valid
        if (is.status == intersection::status::e_inside) {
            is.surface = sf;
            is.direction = std::signbit(vector::dot(is.p3, h.dir(s)))
                               ? intersection::direction::e_opposite
                               : intersection::direction::e_along;
            is.volume_link = mask.volume_link();
        }

        return ret;
    }
};

}  // namespace detail

/// Specialization of the @c helix_intersector for 2D cylindrical surfaces
template <typename transform3_t, typename mask_t>
struct helix_intersector<
    transform3_t, mask_t,
    std::enable_if_t<std::is_same_v<typename mask_t::measurement_frame_type,
                                    cylindrical2<transform3_t>>,
                     void>>
    : public detail::helix_cylinder_intersector<transform3_t> {
    using intersector_impl = detail::helix_cylinder_intersector<transform3_t>;

    using scalar_type = typename intersector_impl::scalar_type;
    using matrix_operator = typename intersector_impl::matrix_operator;
    using point3 = typename intersector_impl::point3;
    using vector3 = typename intersector_impl::vector3;
    using helix_type = typename intersector_impl::helix_type;
};

}  // namespace detray