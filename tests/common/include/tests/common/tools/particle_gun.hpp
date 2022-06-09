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
#include "detray/intersection/cylinder_intersector.hpp"
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

    using intersection_type = line_plane_intersection;

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

        std::vector<std::pair<dindex, intersection_type>> intersection_record;

        // Loop over volumes
        for (const auto &volume : detector.volumes()) {
            for (const auto &[sf_idx, sf] :
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
                    if (sf.is_portal()) {
                        // std::cout << "Found portal: " << sf_idx << std::endl;
                    } else {
                        // std::cout << "Found surface: " << sf_idx <<
                        // std::endl;
                    }
                    // surface the candidate belongs to
                    sfi.index = volume.index();
                    intersection_record.emplace_back(sf_idx, sfi);
                }
            }
        }

        // Sort intersections by distance to origin of the trajectory
        auto sort_path = [&](std::pair<dindex, intersection_type> a,
                             std::pair<dindex, intersection_type> b) -> bool {
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

            // Check all masks of this surface for intersection with the helix.
            // In order to not recompile the geometry, use external intersectors
            for (const auto &mask : range(mask_group, mask_range)) {
                intersection_type sfi;
                // Make sure helix_intersector is only called for the correct
                // surface type
                if constexpr (std::is_same_v<
                                  typename mask_defs::template get_type<
                                      id>::type::intersector_type,
                                  planar_intersector>) {
                    sfi = helix_plane_intersector(ctf, h, mask);
                }
                if constexpr (std::is_same_v<
                                  typename mask_defs::template get_type<
                                      id>::type::intersector_type,
                                  cylinder_intersector>) {
                    sfi = helix_cylinder_intersector(ctf, h, mask);
                }
                // Compare with the ray that keeps only intersection along its
                // direction
                if (sfi.status == intersection::status::e_inside and
                    sfi.direction == intersection::direction::e_along) {
                    sfi.index = volume_index;
                    return sfi;
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
        return intersection_type{};
    }

    /// @brief Intersection implementation for helical trajectories with planar
    /// surfaces.
    ///
    /// The algorithm uses the Newton-Raphson method to find an intersection on
    /// the unbounded surface and then applies the mask.
    ///
    /// @return the intersection.
    template <typename transform_t, typename mask_t,
              std::enable_if_t<std::is_same_v<typename mask_t::intersector_type,
                                              planar_intersector>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline static auto helix_plane_intersector(
        const transform_t &trf, const helix &h, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) -> intersection_type {

        using local_frame = typename mask_t::local_type;

        // Get the surface info
        const auto &sm = trf.matrix();
        auto sn = getter::vector<3>(sm, 0, 2);
        auto st = getter::vector<3>(sm, 0, 3);

        // Starting point on the helix for the Newton iteration
        scalar epsilon{1e-4};
        scalar s{getter::norm(sn) - scalar{0.1}};
        scalar s_prev{s - scalar{0.1}};

        // Guard against inifinite loops
        std::size_t n_tries{0};
        std::size_t max_n_tries{100};

        // f(s) = sn * (h.pos(s) - st) == 0
        // Run the iteration on s
        while (std::abs(s - s_prev) > epsilon and n_tries < max_n_tries) {
            // f'(s) = sn * h.dir(s)
            scalar denom{vector::dot(sn, h.dir(s))};
            if (denom == 0.) {
                break;
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
        is.path = getter::norm(h.pos(s));
        is.p3 = h.pos(s);
        constexpr local_frame local_converter{};
        is.p2 = local_converter(trf, is.p3);
        is.status = mask.template is_inside<local_frame>(is.p2, tolerance);
        // std::cout << "DIRECTION: " << vector::dot(st, h.pos(s)) << std::endl;
        is.direction = vector::dot(st, h.dir(s)) > 0.
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();
        return is;
    }

    /// @brief Intersection implementation for helical trajectories with
    /// cylinder surfaces.
    ///
    /// The algorithm uses the Newton-Raphson method to find an intersection on
    /// the unbounded surface and then applies the mask.
    ///
    /// @return the intersection.
    template <typename transform_t, typename mask_t,
              std::enable_if_t<std::is_same_v<typename mask_t::intersector_type,
                                              cylinder_intersector>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline static auto helix_cylinder_intersector(
        const transform_t &trf, const helix &h, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) -> intersection_type {

        // using point3 = __plugin::point3<scalar>;
        using vector3 = __plugin::vector3<scalar>;
        using local_frame = typename mask_t::local_type;

        // Get the surface info
        const auto &sm = trf.matrix();
        // Cylinder z axis
        vector3 sz = getter::vector<3>(sm, 0, 2);
        // std::cout << "Cylinder axis: " << sz[0] << ", " << sz[1] << ", " <<
        // sz[2] <<  std::endl;

        // Cylinder centre
        auto sc = getter::vector<3>(sm, 0, 3);
        // std::cout << "Cylinder center: " << sc[0] << ", " << sc[1] << ", " <<
        // sc[2] <<  std::endl;

        // Starting point on the helix for the Newton iteration
        scalar epsilon{1e-4};
        // The mask is a cylinder type -> it provides a radius
        scalar r{mask[0]};
        scalar s{r * getter::perp(h.dir(epsilon))};
        // scalar s{r - scalar{0.1}};
        scalar s_prev{s - scalar{0.1}};
        // std::cout << "Start s: " << s << std::endl;
        // std::cout << "Helix start pos: " << h.pos(s)[0] << ", " <<
        // h.pos(s)[1] << ", " << h.pos(s)[2] <<  std::endl; std::cout << "Helix
        // start dir: " << h.dir(s)[0] << ", " << h.dir(s)[1] << ", " <<
        // h.dir(s)[2] <<  std::endl;

        // Guard against inifinite loops
        std::size_t n_tries{0};
        std::size_t max_n_tries{1000};

        // f(s) = ((h.pos(s) - sc) x sz)^2 - r^2 == 0
        // Run the iteration on s
        while (std::abs(s - s_prev) > epsilon and n_tries < max_n_tries) {
            // f'(s) = 2 * ( (h.pos(s) - sc) x sz) * (h.dir(s) x sz) )
            /*std::cout << "pos x ez: " << vector::cross(h.pos(s) - sc, sz) <<
            std::endl;
            std::cout << "dir x ez: " << vector::cross(h.dir(s), sz) <<
            std::endl;*/
            // scalar denom{scalar{2.} * vector::dot(sr, h.dir(s)) *
            // vector::dot(sr, h.pos(s) - sc)};
            scalar denom{scalar{2} *
                         vector::dot(vector::cross(h.pos(s) - sc, sz),
                                     vector::cross(h.dir(s), sz))};
            // std::cout << "denom: " << denom << std::endl;
            if (denom == 0.) {
                break;
            }
            // x_n+1 = x_n - f(s) / f'(s)
            s_prev = s;
            vector3 crp = vector::cross(h.pos(s) - sc, sz);
            s -= (vector::dot(crp, crp) - r * r) / denom;
            // std::cout << "New s: " << s << std::endl;
            // std::cout << "New pos: " << h.pos(s)[0] << ", " << h.pos(s)[1] <<
            // ", " << h.pos(s)[2] <<  std::endl;
            ++n_tries;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return intersection_type{};
        }

        // Build intersection struct from helix parameter s
        intersection_type is;
        is.path = getter::norm(h.pos(s));
        is.p3 = h.pos(s);
        constexpr local_frame local_converter{};
        is.p2 = local_converter(trf, is.p3);
        // std::cout << "local2: " << is.p2[0] << ", " << is.p2[1] << ", " <<
        // is.p2[2] << std::endl;
        auto local3 = trf.point_to_local(is.p3);
        // std::cout << "local3: " << local3[0] << ", " << local3[1] << ", " <<
        // local3[2] << std::endl;
        is.status = mask.template is_inside<local_frame>(local3, tolerance);
        // Explicitely check for radial match
        scalar is_r = getter::perp(is.p3);
        // std::cout << "r: " << is_r << ", " << mask[0] << std::endl;
        if (std::abs(is_r - mask[0]) >= epsilon) {
            // std::cout << "missed portal: " << is_r << ", " << mask[0] <<
            // std::endl;
            is.status = intersection::status::e_missed;
        }
        // std::cout << std::boolalpha << (is.status ==
        // intersection::status::e_inside) << std::endl; std::cout << "Helix
        // final pos: " << is.p3[0] << ", " << is.p3[1] << ", " << is.p3[2] <<
        // std::endl; std::cout << "Helix final dir: " << h.dir(s)[0] << ", " <<
        // h.dir(s)[1] << ", " << h.dir(s)[2] << std::endl;
        is.direction = vector::dot(is.p3, h.dir(s)) > 0.
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();
        return is;
    }
};

}  // namespace detray