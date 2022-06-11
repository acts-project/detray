/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/propagator/track.hpp"

namespace detray {

/// @brief Intersection implementation for helical trajectories with planar
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
struct helix_plane_intersector {

    using intersection_type = line_plane_intersection;

    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;
    using point2 = __plugin::point2<detray::scalar>;
    using cylindrical2 = __plugin::cylindrical2<detray::scalar>;

    /// TODO: Placeholder
    const vector3 default_b_field{0., 0., 2.};

    /// Intersection method for planar surfaces
    ///
    /// @tparam transform_t The type of placement matrix of the surface
    /// @tparam track_t The type of the track (which carries the context
    ///         object)
    /// @tparam mask_t The mask type applied to the local frame
    /// @tparam local_frame The local frame type to be intersected
    ///
    /// Contextual part:
    /// @param trf the surface to be intersected
    /// @param track the track information
    ///
    /// Non-contextual part:
    /// @param mask the local mask
    /// @param tolerance is the mask specific tolerance
    ///
    /// @return the intersection with optional parameters
    template <typename transform_t, typename mask_t,
              std::enable_if_t<std::is_class_v<typename mask_t::local_type>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline auto intersect(
        const transform_t &trf, const free_track_parameters &track,
        const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon,
        vector3 const *const b_field = nullptr) -> intersection_type const {

        return intersect(trf, detail::helix(track, b_field), mask, tolerance);
    }

    /// Intersection method for cylindrical surfaces using Newton-Raphson method
    ///
    /// @tparam transform_t The type of placement matrix of the surface
    /// @tparam mask_t The mask type applied to the local frame
    ///
    /// Contextual part:
    /// @param trf the transform of the surface to be intersected
    ///
    /// Non-contextual part:
    /// @param mask the local mask
    /// @param tolerance is the mask specific tolerance
    ///
    /// @return the intersection with optional parameters
    template <typename transform_t, typename mask_t,
              std::enable_if_t<std::is_class_v<typename mask_t::local_type>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline static auto intersect(
        const transform_t &trf, const detail::helix &h, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) -> intersection_type {

        using local_frame = typename mask_t::local_type;

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
                return intersection_type{};
            }
            // x_n+1 = x_n - f(s) / f'(s)
            s_prev = s;
            s -= vector::dot(sn, h.pos(s) - st) / denom;
            ++n_tries;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return intersection_type{};
        }

        // Build intersection struct from helix parameter s
        intersection_type is;
        const point3 helix_pos = h.pos(s);

        is.path = getter::norm(helix_pos);
        is.p3 = helix_pos;
        constexpr local_frame local_converter{};
        is.p2 = local_converter(trf, is.p3);

        is.status = mask.template is_inside<local_frame>(is.p2, tolerance);
        is.direction = vector::dot(st, h.dir(s)) > scalar{0.}
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();

        return is;
    }
};

}  // namespace detray
