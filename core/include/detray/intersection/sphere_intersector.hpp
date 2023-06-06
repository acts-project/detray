/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <iostream>
#include <type_traits>

namespace detray {

/// A functor to find intersections between straight line and planar surface
template <typename intersection_t>
struct sphere_intersector {

    /// linear algebra types
    /// @{
    using transform3_type = typename intersection_t::transform3_type;
    using scalar_type = typename transform3_type::scalar_type;
    using value_type = typename transform3_type::value_type;
    using point3 = typename transform3_type::point3;
    using point2 = typename transform3_type::point2;
    using vector3 = typename transform3_type::vector3;
    /// @}

    using intersection_type = intersection_t;
    using ray_type = detray::ray<transform3_type>;
#if(IS_SOA)
    /// Operator function to find intersections between ray and planar mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    ///
    /// @return the intersection
    template <typename mask_t, typename surface_t>
    DETRAY_HOST_DEVICE inline std::array<intersection_t, 2> operator()(
        const ray_type &ray, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        intersection_t is;

        const auto r{mask[mask_t::shape::e_r]};
        const vector3 center = trf.translation();

        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();

        const point3 oc = ro - center;
        auto a{vector::dot(rd, rd)};
        auto b{2.f * vector::dot(oc, rd)};
        auto c{vector::dot(oc, oc) - (r * r)};

        const auto qe =
            detail::quadratic_equation<decltype(a)>{a, b, c, 0.f};

        std::array<intersection_t, 2> ret;
        //switch (qe.solutions()) {
        //    case 2:
                ret[1] = build_candidate(ray, mask, trf, qe.larger());
                ret[1].surface = sf;
                // If there are two solutions, reuse the case for a single
                // solution to setup the intersection with the smaller path
                // in ret[0]
        //        [[fallthrough]];
        //    case 1:
                ret[0] = build_candidate(ray, mask, trf, qe.smaller());
                ret[0].surface = sf;
        //        break;
        //    case 0:
        //        ret[0].status = intersection::status::e_missed;
        //        ret[1].status = intersection::status::e_missed;
        //};

        // Even if there are two geometrically valid solutions, the smaller one
        // might not be passed on if it is below the overstepping tolerance:
        // see 'build_candidate'
        return ret;
    }

    /// From the intersection path, construct an intersection candidate and
    /// check it against the surface boundaries (mask).
    ///
    /// @returns the intersection candidate. Might be (partially) uninitialized
    /// if the overstepping tolerance is not met or the intersection lies
    /// outside of the mask.
    template <typename mask_t, typename scalar_t>
    DETRAY_HOST_DEVICE inline intersection_t build_candidate(
        const ray_type &ray, const mask_t &mask, const transform3_type &trf,
        const scalar_t path) const {

        intersection_t is;

        // Construct the candidate only when needed
        //if (path >= ray.overstep_tolerance()) {

            const point3 &ro = ray.pos();
            const vector3 &rd = ray.dir();

            is.path = path;
            const point3 p3 = ro + is.path * rd;

            // No further mask check needed, if the quadratic equation found a
            // solution, an intersection is guaranteed
            /*is.local = mask.to_local_frame(trf, p3);
            is.status = intersection::status::e_inside;

            is.direction = std::signbit(is.path)
                               ? intersection::direction::e_opposite
                               : intersection::direction::e_along;
            is.volume_link = mask.volume_link();

            // Get incidence angle
            const vector3 normal = mask.normal(is.local);
            is.cos_incidence_angle = vector::dot(rd, normal);*/
        //} else {
        //    is.status = intersection::status::e_missed;
        //}

        return is;
    }
#else
    /// Operator function to find intersections between ray and planar mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    ///
    /// @return the intersection
    template <typename mask_t, typename surface_t>
    DETRAY_HOST_DEVICE inline std::array<intersection_t, 2> operator()(
        const ray_type &ray, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        intersection_t is;

        const auto r{mask[mask_t::shape::e_r]};
        const vector3 center = trf.translation();

        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();

        const point3 oc = ro - center;
        auto a{vector::dot(rd, rd)};
        auto b{2.f * vector::dot(oc, rd)};
        auto c{vector::dot(oc, oc) - (r * r)};

        const auto qe =
            detail::quadratic_equation<decltype(a)>{a, b, c, 0.f};

        std::array<intersection_t, 2> ret;
        switch (qe.solutions()) {
            case 2:
                ret[1] = build_candidate(ray, mask, trf, qe.larger());
                ret[1].surface = sf;
                // If there are two solutions, reuse the case for a single
                // solution to setup the intersection with the smaller path
                // in ret[0]
                [[fallthrough]];
            case 1:
                ret[0] = build_candidate(ray, mask, trf, qe.smaller());
                ret[0].surface = sf;
                break;
            case 0:
                ret[0].status = intersection::status::e_missed;
                ret[1].status = intersection::status::e_missed;
        };

        // Even if there are two geometrically valid solutions, the smaller one
        // might not be passed on if it is below the overstepping tolerance:
        // see 'build_candidate'
        return ret;
    }

    /// From the intersection path, construct an intersection candidate and
    /// check it against the surface boundaries (mask).
    ///
    /// @returns the intersection candidate. Might be (partially) uninitialized
    /// if the overstepping tolerance is not met or the intersection lies
    /// outside of the mask.
    template <typename mask_t, typename scalar_t>
    DETRAY_HOST_DEVICE inline intersection_t build_candidate(
        const ray_type &ray, const mask_t &mask, const transform3_type &trf,
        const scalar_t path) const {

        intersection_t is;

        // Construct the candidate only when needed
        if (path >= ray.overstep_tolerance()) {

            const point3 &ro = ray.pos();
            const vector3 &rd = ray.dir();

            is.path = path;
            const point3 p3 = ro + is.path * rd;

            // No further mask check needed, if the quadratic equation found a
            // solution, an intersection is guaranteed
            is.local = mask.to_local_frame(trf, p3);
            is.status = intersection::status::e_inside;

            is.direction = std::signbit(is.path)
                               ? intersection::direction::e_opposite
                               : intersection::direction::e_along;
            is.volume_link = mask.volume_link();

            // Get incidence angle
            const vector3 normal = mask.normal(is.local);
            is.cos_incidence_angle = vector::dot(rd, normal);
        } else {
            is.status = intersection::status::e_missed;
        }

        return is;
    }
#endif
    /// Operator function to updtae an intersections between a ray and a planar
    /// surface.
    ///
    /// @tparam mask_t is the input mask type
    ///
    /// @param ray is the input ray trajectory
    /// @param sfi the intersection to be updated
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray, intersection_t &sfi, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {
        sfi = this->operator()(ray, sfi.surface, mask, trf, mask_tolerance);
    }
};

}  // namespace detray
