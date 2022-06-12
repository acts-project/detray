/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>

#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"

namespace detray {

/// @brief Intersection implementation for ray trajectories with planar
/// surfaces.
struct ray_plane_intersector {

    using intersection_type = line_plane_intersection;

    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /// Intersection method for planar surfaces
    ///
    /// @tparam transform_t The type of placement matrix of the cylinder surface
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
    template <
        typename transform_t, typename track_t, typename mask_t,
        std::enable_if_t<std::is_class_v<typename mask_t::local_type>, bool> =
            true,
        std::enable_if_t<not std::is_same_v<track_t, detail::ray>, bool> = true>
    DETRAY_HOST_DEVICE inline intersection_type intersect(
        const transform_t &trf, const track_t &track, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) const {

        return intersect(trf, detail::ray(track), mask, tolerance,
                         track.overstep_tolerance());
    }

    /// Intersection method for planar surfaces
    ///
    /// @tparam transform_t The type of placement matrix of the cylinder surface
    /// @tparam mask_t The mask type applied to the local frame
    /// @tparam local_frame The local frame type to be intersected
    ///
    /// Contextual part:
    /// @param trf the transform of the surface to be intersected
    /// @param ro the origin of the ray
    /// @param rd the direction of the ray
    ///
    /// Non-contextual part:
    /// @param mask the local mask
    /// @param tolerance is the mask specific tolerance
    ///
    /// @return the intersection with optional parameters
    template <typename transform_t, typename mask_t,
              std::enable_if_t<std::is_class_v<typename mask_t::local_type>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline intersection_type intersect(
        const transform_t &trf, const detail::ray ray, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon,
        const scalar overstep_tolerance = 0.) const {

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
            intersection_type is;
            is.path = vector::dot(sn, st - ro) / denom;
            is.p3 = ro + is.path * rd;
            constexpr local_frame local_converter{};
            is.p2 = local_converter(trf, is.p3);
            is.status = mask.template is_inside<local_frame>(is.p2, tolerance);
            is.direction = is.path > overstep_tolerance
                               ? intersection::direction::e_along
                               : intersection::direction::e_opposite;
            is.link = mask.volume_link();
            return is;
        }
        return intersection_type{};
    }
};

}  // namespace detray
