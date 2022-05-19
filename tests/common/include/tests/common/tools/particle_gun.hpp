/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// system include
#include <cmath>
#include <iostream>
#include <utility>

// detray include(s)
#include "detray/definitions/units.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/intersection/planar_intersector.hpp"
#include "detray/utils/enumerate.hpp"
#include "tests/common/tools/test_trajectories.hpp"

namespace detray {

/// @brief struct that holds functionality to shoot a ray through a detector.
///
/// Records intersections with every detector surface along the ray.
struct particle_gun {

    using intersection_t = line_plane_intersection;

    /// Intersect all surfaces in a detector with a given ray.
    ///
    /// @param detector the detector.
    /// @param traj the trajectory to be shot through the detector.
    ///
    /// @return a sorted vector of volume indices with the corresponding
    ///         intersections of the surfaces that were encountered.
    template <typename detector_t, typename trajectory_t>
    DETRAY_HOST_DEVICE inline static auto shoot_particle(
        const detector_t &detector, const trajectory_t &traj) {

        std::vector<std::pair<dindex, intersection_t>> intersection_record;

        // Loop over volumes
        for (const auto &volume : detector.volumes()) {
            for (const auto [sf_idx, sf] :
                 enumerate(detector.surfaces(), volume)) {
                // Retrieve candidate from the surface
                auto sfi = intersect(traj, sf, detector.transform_store(),
                                     detector.mask_store());

                // Candidate is invalid if it oversteps too far (this is neg!)
                if (sfi.path < traj.overstep_tolerance()) {
                    continue;
                }
                // Accept if inside
                if (sfi.status == intersection::status::e_inside) {
                    // surface the candidate belongs to
                    sfi.index = volume.index();
                    intersection_record.emplace_back(sf_idx, sfi);
                }
            }
        }

        // Sort intersections by distance to origin of the trajectory
        auto sort_path = [&](std::pair<dindex, intersection_t> a,
                             std::pair<dindex, intersection_t> b) -> bool {
            return (a.second < b.second);
        };
        std::sort(intersection_record.begin(), intersection_record.end(),
                  sort_path);

        return intersection_record;
    }

    /// Wrap the @c intersection_kernel intersect function for rays
    template <typename surface_t, typename transform_container,
              typename mask_container>
    DETRAY_HOST_DEVICE static inline auto intersect(
        const ray &r, surface_t &surface,
        const transform_container &contextual_transforms,
        const mask_container &masks) {
        return detray::intersect(r, surface, contextual_transforms, masks);
    }

    /// Start helix intersection unrolling. See documentation of the detray
    /// intersection kernel
    template <typename surface_t, typename transform_container,
              typename mask_container>
    DETRAY_HOST_DEVICE static inline auto intersect(
        const helix &h, surface_t &surface,
        const transform_container &contextual_transforms,
        const mask_container &masks) {
        // Gather all information to perform intersections
        const auto &ctf = contextual_transforms[surface.transform()];
        const auto volume_index = surface.volume();
        const auto mask_id = surface.mask_type();
        const auto &mask_range = surface.mask_range();

        // No intersectors available for non-planar surfaces
        if (surface.is_portal()) {
            return intersection_t{};
        }
        // Unroll the intersection depending on the mask container size
        using mask_defs = typename surface_t::mask_defs;

        return unroll_helix_intersect<mask_defs>(
            h, ctf, masks, mask_range, mask_id, volume_index,
            std::make_integer_sequence<unsigned int, mask_defs::n_types>{});
    }

    /// Helix version of the intersection unrolling. Calls the
    /// @c helix_intersector of this class instead of the mask's native
    /// intersector.
    template <typename mask_defs, typename transform_t,
              typename mask_container_t, typename mask_range_t,
              unsigned int first_mask_id, unsigned int... remaining_mask_ids>
    DETRAY_HOST_DEVICE static inline auto unroll_helix_intersect(
        const helix &h, const transform_t &ctf, const mask_container_t &masks,
        const mask_range_t &mask_range,
        const typename mask_container_t::id_type mask_id, dindex volume_index,
        std::integer_sequence<unsigned int, first_mask_id,
                              remaining_mask_ids...>
        /*available_ids*/) {

        // Pick the first one for interseciton
        if (mask_id == first_mask_id) {
            // Get the mask id that was found
            constexpr auto id = mask_container_t::to_id(first_mask_id);
            auto &mask_group = masks.template group<id>();

            // Check all masks of this surface for intersection with the helix
            for (const auto &mask : range(mask_group, mask_range)) {
                // Make sure helix_intersector is only called for planar surface
                if constexpr (std::is_same_v<
                                  typename mask_defs::template get_type<
                                      id>::type::intersector_type,
                                  planar_intersector>) {
                    std::cout << "Intersecting planar sf with helix"
                              << std::endl;
                    auto sfi = helix_intersector(ctf, h, mask);
                    std::cout << sfi.p3[0] << ", " << sfi.p3[1] << ", "
                              << sfi.p3[2] << std::endl;
                    if (sfi.status == intersection::status::e_inside) {
                        sfi.index = volume_index;
                        return sfi;
                    }
                }
            }
        }

        // The reduced integer sequence
        std::integer_sequence<unsigned int, remaining_mask_ids...> remaining;

        // Unroll as long as you have at least 1 entries
        if constexpr (remaining.size() >= 1) {
            return (unroll_helix_intersect<mask_defs>(
                h, ctf, masks, mask_range, mask_id, volume_index, remaining));
        }

        // No intersection was found
        return intersection_t{};
    }

    /// @brief Intersection implementation for helical trajectories.
    ///
    /// The algorithm uses the Newton-Raphson method to find an intersection on
    /// the unbounded surface and then applies the mask.
    ///
    /// @return the intersection.
    template <typename transform_t, typename mask_t,
              std::enable_if_t<std::is_same_v<typename mask_t::intersector_type,
                                              planar_intersector>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline static auto helix_intersector(
        const transform_t &trf, const helix &h, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) -> intersection_t {

        using local_frame = typename mask_t::local_type;

        // Get the surface info
        const auto &sm = trf.matrix();
        auto sn = getter::vector<3>(sm, 0, 2);
        auto st = getter::vector<3>(sm, 0, 3);

        // starting point on the helix for the Newton iteration
        scalar epsilon = 1e-7;
        scalar s{getter::norm(sn) - 0.1};
        scalar s_prev{s - 0.2};
        std::cout << "Starting point: " << s << std::endl;

        // Guard against inifinite loops
        std::size_t n_tries{0};
        std::size_t max_n_tries = 1000;
        const auto dir = h.pos(9.);
        // std::cout << "normal: " << sn[0] << ", " << sn[1] << ", " << sn[2] <<
        // std::endl; std::cout << "dir: " << dir[0] << ", " << dir[1] << ", "
        // << dir[2] << std::endl;
        //  Run the iteration on s
        while (std::abs(s - s_prev) > epsilon and n_tries < max_n_tries) {
            scalar denom = vector::dot(sn, h.dir(s));
            if (denom == 0.) {
                std::cout << "denom zero!" << std::endl;
                break;
            }
            s_prev = s;
            // std::cout << "denom: " << denom << std::endl;
            s -= vector::dot(sn, h.pos(s) - st) / denom;
            ++n_tries;

            std::cout << "Step: " << n_tries << ", s: " << s << std::endl;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return intersection_t{};
        }

        // Build intersection struct from helix parameter s
        intersection_t is;
        is.path = getter::norm(h.pos(s));
        is.p3 = h.pos(s);
        constexpr local_frame local_converter{};
        is.p2 = local_converter(trf, is.p3);
        is.status = mask.template is_inside<local_frame>(is.p2, tolerance);
        is.direction = is.path > h.overstep_tolerance()
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();
        return is;
    }
};

}  // namespace detray
