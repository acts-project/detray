/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical2D.hpp"
#include "detray/definitions/boolean.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/soa/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray::soa {

/// A functor to find intersections between a ray and a 2D cylinder mask
template <typename intersection_t>
struct cylinder_intersector {

    /// Linear algebra types
    /// @{
    using algebra = typename intersection_t::algebra;
    using transform3_type = typename intersection_t::transform3D;
    using value_type = typename intersection_t::value_t;
    using scalar_type = typename intersection_t::scalar_t;
    using point3 = typename intersection_t::point3D;
    using point2 = typename intersection_t::point2D;
    using vector3 = typename intersection_t::vector3D;
    /// @}

    using intersection_type = intersection_t;
    using ray_type =
        detray::detail::ray<dtransform3D<ALGEBRA_PLUGIN<value_type>>>;

    /// Operator function to find intersections between a ray and a 2D cylinder
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
    /// @return the intersections.
    template <typename mask_t, typename surface_t,
              std::enable_if_t<std::is_same_v<typename mask_t::local_frame_type,
                                              cylindrical2D<algebra>>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline std::array<intersection_t, 2> operator()(
        const ray_type &ray, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        // One or both of these solutions might be invalid
        const auto qe = solve_intersection(ray, mask, trf);

        std::array<intersection_t, 2> ret;
        ret[1] = build_candidate(ray, mask, trf, qe.larger(), mask_tolerance);
        ret[1].sf_desc = sf;

        ret[0] = build_candidate(ray, mask, trf, qe.smaller(), mask_tolerance);
        ret[0].sf_desc = sf;

        // Even if there are two geometrically valid solutions, the smaller one
        // might not be passed on if it is below the overstepping tolerance:
        // see 'build_candidate'
        return ret;
    }

    /// Operator function to find intersections between a ray and a 2D cylinder
    ///
    /// @tparam mask_t is the input mask type
    ///
    /// @param ray is the input ray trajectory
    /// @param sfi the intersection to be updated
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    template <typename mask_t,
              std::enable_if_t<
                  std::is_same_v<typename mask_t::local_frame_type,
                                 cylindrical2D<ALGEBRA_PLUGIN<scalar_type>>>,
                  bool> = true>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray, intersection_t &sfi, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        // One or both of these solutions might be invalid
        const auto qe = solve_intersection(ray, mask, trf);

        // Construct the candidate only when needed
        sfi.status = (qe.solutions() > 0.f);

        if (detray::detail::none_of(sfi.status)) {
            return sfi;
        }

        sfi = build_candidate(ray, mask, trf, qe.smaller(), mask_tolerance);
    }

    protected:
    /// Calculates the distance to the (two) intersection points on the
    /// cylinder in global coordinates.
    ///
    /// @returns a quadratic equation object that contains the solution(s).
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline auto solve_intersection(
        const ray_type &ray, const mask_t &mask,
        const transform3_type &trf) const {
        const auto &m = trf.matrix();
        const vector3 sz = getter::vector<3>(m, 0u, 2u);
        const vector3 sc = getter::vector<3>(m, 0u, 3u);

        const scalar_type r = mask[mask_t::shape::e_r];

        const auto &pos = ray.pos();
        const auto &dir = ray.dir();
        const point3 ro{pos[0], pos[1], pos[2]};
        const vector3 rd{dir[0], dir[1], dir[2]};

        const vector3 tmp = ro - sc;
        const auto pc_cross_sz = vector::cross(tmp, sz);
        const auto rd_cross_sz = vector::cross(rd, sz);
        const scalar_type a = vector::dot(rd_cross_sz, rd_cross_sz);
        const scalar_type b = 2.f * vector::dot(rd_cross_sz, pc_cross_sz);
        const scalar_type c = vector::dot(pc_cross_sz, pc_cross_sz) - (r * r);

        return soa::detail::quadratic_equation{a, b, c};
    }

    /// From the intersection path, construct an intersection candidate and
    /// check it against the surface boundaries (mask).
    ///
    /// @returns the intersection candidate. Might be (partially) uninitialized
    /// if the overstepping tolerance is not met or the intersection lies
    /// outside of the mask.
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_t build_candidate(
        const ray_type &ray, const mask_t &mask, const transform3_type &trf,
        const scalar_type path, const scalar_type mask_tolerance = 0.f) const {

        intersection_t is;

        const auto &pos = ray.pos();
        const auto &dir = ray.dir();
        const point3 ro{pos[0], pos[1], pos[2]};
        const vector3 rd{dir[0], dir[1], dir[2]};

        is.path = path;
        const point3 p3 = ro + is.path * rd;

        is.local = mask.to_local_frame(trf, p3);
        is.status = mask.is_inside(is.local, mask_tolerance);

        is.direction = math_ns::signbit(is.path);
        is.volume_link = mask.volume_link();

        // Get incidence angle
        const scalar_type phi{is.local[0] / is.local[2]};
        const vector3 normal = {math_ns::cos(phi), math_ns::sin(phi), 0.f};
        is.cos_incidence_angle = vector::dot(rd, normal);

        // Mask the values where the overstepping tolerance was not met
        is.status &= (is.path >= ray.overstep_tolerance());

        return is;
    }
};

}  // namespace detray::soa
