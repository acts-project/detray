/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>

#include "detray/core/intersection.hpp"
#include "detray/core/track.hpp"

namespace detray {

/** This is an intersector struct for planar surface
 */
struct planar_intersector {

    using transform3 = __plugin::transform3;
    using point3 = __plugin::point3;
    using vector3 = __plugin::vector3;
    using point2 = __plugin::point2;

    /** Intersection method for planar surfaces
     *
     * @tparam track_type The type of the track (which carries the context
     *         object)
     * @tparam mask_type The mask type applied to the local frame
     * @tparam local_frame The local frame type to be intersected
     *
     * Contextual part:
     * @param trf the surface to be intersected
     * @param track the track information
     *
     * Non-contextual part:
     * @param mask the local mask
     * @param tolerance is the mask specific tolerance
     *
     * @return the intersection with optional parameters
     **/
    template <typename track_type, typename mask_type,
              std::enable_if_t<std::is_class_v<typename mask_type::local_type>,
                               bool> = true>
    inline intersection intersect(
        const transform3 &trf, const track_type &track, const mask_type &mask,
        const typename mask_type::mask_tolerance &tolerance =
            mask_type::within_epsilon) const {
        return intersect(trf, track.pos, track.dir, mask, tolerance,
                         track.overstep_tolerance);
    }

    /** Intersection method for planar surfaces
     *
     * @tparam mask_type The mask type applied to the local frame
     * @tparam local_frame The local frame type to be intersected
     *
     * Contextual part:
     * @param trf the transform of the surface to be intersected
     * @param ro the origin of the ray
     * @param rd the direction of the ray
     *
     * Non-contextual part:
     * @param mask the local mask
     * @param tolerance is the mask specific tolerance
     *
     * @return the intersection with optional parameters
     **/
    template <typename mask_type,
              std::enable_if_t<std::is_class_v<typename mask_type::local_type>,
                               bool> = true>
    inline intersection intersect(const transform3 &trf, const point3 &ro,
                                  const vector3 &rd, const mask_type &mask,
                                  const typename mask_type::mask_tolerance
                                      &tolerance = mask_type::within_epsilon,
                                  scalar overstep_tolerance = 0.) const {

        using local_frame = typename mask_type::local_type;

        // Retrieve the surface normal & translation (context resolved)
        const auto &sm = trf.matrix();
        auto sn = getter::vector<3>(sm, 0, 2);
        auto st = getter::vector<3>(sm, 0, 3);

        // Intersection code
        scalar denom = vector::dot(rd, sn);
        if (denom != 0.0) {
            intersection is;
            is.path = vector::dot(sn, st - ro) / denom;
            is.p3 = ro + is.path * rd;
            constexpr local_frame local_converter{};
            is.p2 = local_converter(trf, is.p3);
            is.status = mask.template is_inside<local_frame>(is.p2, tolerance);
            is.direction = is.path > overstep_tolerance ? e_along : e_opposite;
            return is;
        }
        return intersection{};
    }
};

}  // namespace detray
