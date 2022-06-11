/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <cmath>
#include <type_traits>
#include <utility>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/propagator/track.hpp"

namespace detray {

namespace detail {

struct unbound;

}

/// @brief Intersection implementation for cylinder surfaces using helical
/// trajectories.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask. On the @c cylinder3
/// mask, it switches on the check of the radial distance.
struct helix_cylinder_intersector {

    using intersection_type = line_plane_intersection;

    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;
    using cylindrical2 = __plugin::cylindrical2<detray::scalar>;

    const vector3 dummyBfield{0., 0., 2.};

    /// Intersection method for a track with a cylindrical surfaces.
    ///
    /// It biulds a helix trajectory from the track parameters and then dele-
    /// gates to a dedicated helix-cylinder intersection implementation.
    ///
    /// @tparam transform_t The type of placement matrix of the cylinder surface
    /// @tparam mask_t The mask type applied to the local frame
    ///
    /// Contextual part:
    /// @param trf the transform of the surface to be intersected
    /// @param track the track information used in the propagation
    ///
    /// Non-contextual part:
    /// @param mask the local cylindrical mask
    /// @param tolerance is the mask specific tolerance
    ///
    /// @return the intersection with optional parameters
    template <
        typename transform_t, typename mask_t,
        std::enable_if_t<
            std::is_same_v<typename mask_t::local_type, cylindrical2> or
                std::is_same_v<typename mask_t::local_type, detail::unbound>,
            bool> = true>
    DETRAY_HOST_DEVICE inline auto intersect(
        const transform_t &trf, const free_track_parameters &track,
        const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) -> intersection_type const {
        return intersect(trf, detail::helix(track, &dummyBfield), mask,
                         tolerance);
    }

    /// Intersection of a helical trajectory with a cylindrical surface.
    ///
    /// @tparam transform_t The type of placement matrix of the cylinder surface
    /// @tparam mask_t The mask type applied to the local frame
    ///
    /// Contextual part:
    /// @param trf the transform of the surface to be intersected
    /// @param h   the helical trajectory, parametrized by its path length
    ///
    /// Non-contextual part:
    /// @param mask the local cylindrical mask
    /// @param tolerance is the mask specific tolerance
    ///
    /// @return the intersection with optional parameters
    template <
        typename transform_t, typename mask_t,
        std::enable_if_t<
            std::is_same_v<typename mask_t::local_type, cylindrical2> or
                std::is_same_v<typename mask_t::local_type, detail::unbound>,
            bool> = true>
    DETRAY_HOST_DEVICE inline static auto intersect(
        const transform_t &trf, const detail::helix &h, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) -> intersection_type const {

        using local_frame = typename mask_t::local_type;

        // Guard against inifinite loops
        const std::size_t max_n_tries{100};
        // Tolerance for convergence
        const scalar tol{1e-3};

        // Get the surface placement
        const auto &sm = trf.matrix();
        // Cylinder z axis
        const vector3 sz = getter::vector<3>(sm, 0, 2);
        // Cylinder centre
        const point3 sc = getter::vector<3>(sm, 0, 3);

        // Starting point on the helix for the Newton iteration
        // The mask is a cylinder -> it provides its radius as the first value
        scalar r{mask[0]};
        // Helix path length parameter
        scalar s{r * getter::perp(h.dir(tol))};
        // Path length in the previous iteration step
        scalar s_prev{s - scalar{0.1}};

        // f(s) = ((h.pos(s) - sc) x sz)^2 - r^2 == 0
        // Run the iteration on s
        std::size_t n_tries{0};
        while (std::abs(s - s_prev) > tol and n_tries < max_n_tries) {

            // f'(s) = 2 * ( (h.pos(s) - sc) x sz) * (h.dir(s) x sz) )
            const vector3 crp = vector::cross(h.pos(s) - sc, sz);
            const scalar denom{scalar{2.} *
                               vector::dot(crp, vector::cross(h.dir(s), sz))};
            // No intersection can be found if dividing by zero
            if (denom == scalar{0.}) {
                return intersection_type{};
            }
            // x_n+1 = x_n - f(s) / f'(s)
            s_prev = s;
            s -= (vector::dot(crp, crp) - r * r) / denom;

            ++n_tries;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return intersection_type{};
        }

        // Build intersection struct from helix parameter s
        intersection_type is;
        point3 helix_pos = h.pos(s);

        is.path = getter::norm(helix_pos);
        is.p3 = std::move(helix_pos);
        constexpr local_frame local_converter{};
        is.p2 = local_converter(trf, is.p3);

        auto local3 = trf.point_to_local(is.p3);
        // Explicitly check for radial match
        is.status =
            mask.template is_inside<local_frame, true>(local3, tolerance);
        is.direction = vector::dot(is.p3, h.dir(s)) > scalar{0.}
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();

        return is;
    }
};

}  // namespace detray
