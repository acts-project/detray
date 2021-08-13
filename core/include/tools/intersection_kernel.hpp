/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/track.hpp"
#include "core/intersection.hpp"
#include "utils/indexing.hpp"
#include "utils/enumerate.hpp"

#include <utility>
#include <tuple>
#include <typeinfo>

namespace detray
{

    /// Transform definition
    using transform3 = __plugin::transform3;

    /** Specialized method to update an intersection when the mask group is resolved
     * 
     * @tparam track_type is the type of the track information including context 
     * @tparam mask_group is the type of the split out mask group from all masks
     * @tparam mask_range is the type of the according mask range object
     * 
     * @param track is the track information 
     * @param trf is the surface transform for the intersection
     * @param masks is the associated (and split out) mask group
     * @param range is the range list of masks to be processed
     * 
     * @note since all masks are on the same surface, only maximally one mask solution at a time
     * is possbile. Once this solution is found, it is returned
     * 
     * @return the surface intersection 
     **/
        template <dindex current_idx = 0,
                  typename intersection_type,
                  typename track_type,
                  typename mask_container,
                  typename mask_range>
        void unroll_intersect(intersection_type &intersection,
                              const track_type &track,
                              const transform3 &ctf,
                              const mask_container &masks,
                              const mask_range &range,
                              const dindex &mask_context)
        {
            // Intersect only the correct mask type
            if (mask_context == current_idx)
            {
                const auto &surface_masks = std::get<current_idx>(masks);
                for (size_t i = range[0]; i < range[1]; i++)
                {
                    const auto &mask = surface_masks[i];
                    auto local = mask.local();
                    intersection = mask.intersector().intersect(ctf, track, local, mask);
                    if (intersection.status == e_inside) {
                        intersection.index[2] = i;
                        return;
                    }
                }
            }
            // Next mask type
            if constexpr (current_idx < std::tuple_size_v<mask_container> - 1) {
                unroll_intersect<current_idx + 1>(intersection, track, ctf, masks, range, mask_context);
            }
        }

    /** Actual kernel method that updates the updates the intersections
     *
     * @tparam surface_container is The type of the surface container
     * @tparam transform_container is the type of the transform container
     * @tparam mask_container is the type of the type of the mask container
     * @tparam range_type is how ranges in the surface container are encoded
     *
     * @param track the track information including the contexts
     * @param surfaces the surface type to be intersected
     * @param contextual_transforms the transform container
     * @param masks the tuple mask container to for the intersection
     * 
     * @return an intersection and the link to the result surface, relative
     *         to the first surface in all batches
     */
    template <typename surface_container,
              typename transform_container,
              typename mask_container,
              typename range_type>
    const auto
    intersect(const track<typename transform_container::context> &track,
              const surface_container &surfaces,
              const range_type &surface_range,
              const transform_container &contextual_transforms,
              const mask_container &masks) -> std::vector<intersection>
    {
        // Output intersections
        std::vector<intersection> result_intersections;
        result_intersections.reserve((surface_range[1] - surface_range[0]) * surfaces[surface_range[0]].n_surfaces);
        // Run over all surface batches in range (each has potentially a unique
        // mask type)
        for (size_t sfbi = surface_range[0]; sfbi < surface_range[1]; sfbi++)
        {
            const auto &surface_batch = surfaces[sfbi];
            dindex mask_type = surface_batch.mask_type;
            for (size_t si = 0; si < surface_batch.n_surfaces; si++)
            {
                const auto ctf = contextual_transforms[surface_batch.transform_idx + si];
                std::array<dindex, 2> mask_range = surface_batch.mask_range_by_surface(si);
                intersection sfi;
                sfi.index = {sfbi, si, dindex_invalid};

                unroll_intersect(sfi, track, ctf, masks, mask_range, mask_type);
                result_intersections.emplace_back(sfi);
            }
        }
        return result_intersections;
    }

    /** Actual kernel method that updates the updates the intersections
     *
     * @tparam surface_container is The type of the surface container
     * @tparam transform_container is the type of the transform container
     * @tparam mask_container is the type of the type of the mask container
     * @tparam range_type is how ranges in the surface container are encoded
     *
     * @param track the track information including the contexts
     * @param surfaces the surface type to be intersected
     * @param contextual_transforms the transform container
     * @param masks the tuple mask container to for the intersection
     * 
     * @return an intersection and the link to the result surface, relative
     *         to the first surface in all batches
     */
    template <typename surface_container,
              typename transform_container,
              typename mask_container>
    const auto
    intersect(const track<typename transform_container::context> &track,
              const std::array<dindex, 3> &sf_index,
              const surface_container &surfaces,
              const transform_container &contextual_transforms,
              const mask_container &masks)
    {
        const auto &surface_batch = surfaces[sf_index[0]];
        const auto mask_type          = surface_batch.mask_type;
        const auto &ctf = contextual_transforms[surface_batch.transform_idx
                                                + sf_index[1]];
        const auto &mask_range    = surface_batch.mask_range_by_surface(sf_index[1]);

        // Create a return intersection and run the variadic unrolling
        intersection sfi;
        sfi.index = sf_index;

        // Unroll the intersection depending on the mask container size
        unroll_intersect(sfi, track, ctf, masks, mask_range, mask_type);

        // Return the (eventually update) intersection and links
        return sfi;
    }

    template <typename range_type,
              typename transform_container,
              typename mask_container>
    const auto
    intersect(const track<typename transform_container::context> &track,
              const dindex sf_index,
              const dindex &mask_type,
              const range_type &mask_range,
              const transform_container &contextual_transforms,
              const mask_container &masks)
    {
        const auto &ctf = contextual_transforms[sf_index];

        // Create a return intersection and run the variadic unrolling
        intersection sfi;
        sfi.index = {dindex_invalid, sf_index, dindex_invalid};

        // Unroll the intersection depending on the mask container size
        unroll_intersect(sfi, track, ctf, masks, mask_range, mask_type);

        // Return the (eventually update) intersection and links
        return sfi;
    }
}
