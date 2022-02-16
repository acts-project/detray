/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <tuple>
#include <utility>

#include "detray/core/intersection.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/tools/track.hpp"
#include "detray/utils/enumerate.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

/// Transform definition
using transform3 = __plugin::transform3<detray::scalar>;

/** Variadic unrolled intersection - any integer sequence
 *
 * @tparam track_t The type of the surface container
 * @tparam mask_container The type of the type of the mask container
 * @tparam mask_range The mask range type
 * @tparam first_mask_id The first mask group id
 *
 * @param tack the track information containing the context
 * @param ctf the contextual transform (context resolved)
 * @param masks the masks container
 * @param rng the range within the mask group to be checked
 * @param mask_id the current mask group id
 * @param available_ids the mask ids to be checked (only needed to set the
 *                      first mask id for the call)
 *
 * @return an intersection struct (invalid if no intersection was found)
 */
template <typename track_t, typename mask_container, typename mask_range,
          unsigned int first_mask_id, unsigned int... remaining_mask_ids>
DETRAY_HOST_DEVICE inline auto unroll_intersect(
    const track_t &track, const transform3 &ctf, const mask_container &masks,
    const mask_range &rng, const unsigned int mask_id, dindex volume_index,
    std::integer_sequence<unsigned int, first_mask_id, remaining_mask_ids...>
    /*available_ids*/) {

    // Pick the first one for interseciton
    if (mask_id == first_mask_id) {

        auto &mask_group = masks.template group<first_mask_id>();

        // Check all masks of this surface for intersection
        for (const auto &mask : range(mask_group, rng)) {
            auto sfi =
                std::move(mask.intersector().intersect(ctf, track, mask));

            if (sfi.status == e_inside) {
                sfi.index = volume_index;

                // Link to next volume is in first position
                sfi.link = detail::get<0>(mask.links());
                return sfi;
            }
        }
    }

    // The reduced integer sequence
    std::integer_sequence<unsigned int, remaining_mask_ids...> remaining;

    // Unroll as long as you have at least 1 entries
    if constexpr (remaining.size() >= 1) {
        return (unroll_intersect(track, ctf, masks, rng, mask_id, volume_index,
                                 remaining));
    }

    // No intersection was found
    return intersection{};
}

/** Kernel method that updates the intersections
 *
 * @tparam track_t The type of the track/context
 * @tparam surface_t The type of the surface container
 * @tparam transform_container The type of the transform container
 * @tparam mask_container The type of the type of the mask container
 *
 * @param track the track information including the contexts
 * @param surface the surface type to be intersected
 * @param contextual_transform the transform container
 * @param masks the tuple mask container to for the intersection
 *
 * @return  an intersection struct (invalid if no intersection was found)
 **/
template <typename track_t, typename surface_t, typename transform_container,
          typename mask_container>
DETRAY_HOST_DEVICE inline const auto intersect(
    const track_t &track, surface_t &surface,
    const transform_container &contextual_transforms,
    const mask_container &masks) {
    // Gather all information to perform intersections
    const auto &ctf = contextual_transforms[surface.transform()];
    const auto &volume_index = surface.volume();
    const auto &mask_link = surface.mask();
    const auto &mask_id = detail::get<0>(mask_link);
    const auto &mask_range = detail::get<1>(mask_link);

    // Unroll the intersection depending on the mask container size
    return unroll_intersect(
        track, ctf, masks, mask_range, mask_id, volume_index,
        std::make_integer_sequence<
            unsigned int,
            detail::tuple_size<typename mask_container::mask_tuple>::value>{});
}

}  // namespace detray