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
    using transform3_type = typename intersection_t::transform3D;
    using value_type = typename intersection_t::value_t;
    using scalar_type = typename intersection_t::scalar_t;
    using point3 = typename intersection_t::point3D;
    using point2 = typename intersection_t::point2D;
    using vector3 = typename intersection_t::vector3D;
    /// @}

    using intersection_type = intersection_t;
    using ray_type = detray::ray<transform3_type>;
#if (IS_SOA)

    /// Class to solve a quadratic equation of type a * x^2 + b * x + c = 0
    ///
    /// @note If there are no real solutions, the result is undefined
    /// @note The solutions are sorted by default. If there is only one
    /// solution, the larger value is undefined.
    class soa_quadratic_equation {
        public:
        soa_quadratic_equation() = delete;

        /// @brief SoA version
        DETRAY_HOST_DEVICE
        constexpr soa_quadratic_equation(
            const value_type a, const scalar_type &b, const scalar_type &c,
            const scalar_type &tolerance =
                std::numeric_limits<scalar_type>::epsilon()) {

            // Linear case
            if (std::abs(a) <= tolerance[0]) {
                m_solutions = 1.f;
                m_values[0] = -c / b;
            } else {
                const scalar_type discriminant = b * b - (4.f * a) * c;

                const auto two_sol = (discriminant > tolerance);
                const auto one_sol = !two_sol && (discriminant >= 0.f);

                // If there is more than one solution, then a != 0 and q != 0
                if (Vc::any_of(two_sol)) {
                    m_solutions = 2.f;
                    m_solutions.setZeroInverted(two_sol);

                    const scalar_type q =
                        -0.5f * (b + Vc::copysign(Vc::sqrt(discriminant), b));

                    scalar_type first = q / a;
                    scalar_type second = c / q;
                    first.setZeroInverted(two_sol);
                    second.setZeroInverted(two_sol);

                    // Sort the solutions
                    const auto do_swap = Vc::abs(second) <= Vc::abs(first);
                    if (do_swap.isFull()) {
                        m_values = {second, first};
                    } else if (do_swap.isEmpty()) {
                        m_values = {first, second};
                    } else {
                        std::cout << "TODO!!!" << std::endl;
                    }
                }

                // Only one solution and a != 0
                if (Vc::any_of(one_sol)) {
                    scalar_type sol = 1.f;
                    scalar_type result = -0.5f * b / a;
                    sol.setZeroInverted(one_sol);
                    result.setZeroInverted(one_sol);

                    m_solutions += sol;
                    m_values[0] += result;
                }
                // discriminant < 0 is not allowed, since all solutions should
                // be real
            }
        }

        /// Getters for the solution(s)
        /// @{
        constexpr const auto &solutions() const { return m_solutions; }
        constexpr const scalar_type &smaller() const { return m_values[0]; }
        constexpr const scalar_type &larger() const { return m_values[1]; }
        /// @}

        private:
        /// Number of solutions of the equation (needs to be floating point to
        /// apply the masks correctly)
        Vc::float_v m_solutions = 0.f;
        /// The solutions
        std::array<scalar_type, 2> m_values{0.f, 0.f};
    };
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
        const transform3_type &trf, const scalar_type = 0.f) const {

        intersection_t is;

        const scalar_type r = mask[mask_t::shape::e_r];
        const vector3 center = trf.translation();

        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();

        const point3 oc = ro - center;
        const scalar_type a = vector::dot(rd, rd);
        const scalar_type b = 2.f * vector::dot(oc, rd);
        const scalar_type c = vector::dot(oc, oc) - (r * r);

        const auto qe = soa_quadratic_equation{a[0], b, c, 0.f};

        std::array<intersection_t, 2> ret;
        ret[1] = build_candidate(ray, mask, trf, qe.larger(), qe.solutions());
        ret[1].sf_desc = sf;

        ret[0] = build_candidate(ray, mask, trf, qe.smaller(), qe.solutions());
        ret[0].sf_desc = sf;

        return ret;
    }

    /// From the intersection path, construct an intersection candidate and
    /// check it against the surface boundaries (mask).
    ///
    /// @returns the intersection candidate. Might be (partially) uninitialized
    /// if the overstepping tolerance is not met or the intersection lies
    /// outside of the mask.
    template <typename mask_t, typename bool_mask_t>
    DETRAY_HOST_DEVICE inline intersection_t build_candidate(
        const ray_type &ray, const mask_t &mask, const transform3_type &trf,
        const scalar_type path, const bool_mask_t &n_solutions) const {

        intersection_t is;

        // Construct the candidate only when needed
        is.status = (n_solutions > 0.f);

        if (is.status.isEmpty()) {
            return is;
        }

        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();

        is.path = path;
        const point3 p3 = ro + is.path * rd;

        // No further mask check needed, if the quadratic equation found a
        // solution, an intersection is guaranteed
        is.local = mask.to_local_frame(trf, p3);

        is.direction = Vc::isnegative(is.path);
        is.volume_link = mask.volume_link();

        // Get incidence angle
        const vector3 normal = mask.normal(is.local);
        is.cos_incidence_angle = vector::dot(rd, normal);

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

        const auto qe = detail::quadratic_equation<decltype(a)>{a, b, c, 0.f};

        std::array<intersection_t, 2> ret;
        switch (qe.solutions()) {
            case 2:
                ret[1] = build_candidate(ray, mask, trf, qe.larger());
                ret[1].sf_desc = sf;
                // If there are two solutions, reuse the case for a single
                // solution to setup the intersection with the smaller path
                // in ret[0]
                [[fallthrough]];
            case 1:
                ret[0] = build_candidate(ray, mask, trf, qe.smaller());
                ret[0].sf_desc = sf;
                break;
            case 0:
                ret[0].status = false;
                ret[1].status = false;
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
            is.status = true;

            is.direction = !std::signbit(is.path);
            is.volume_link = mask.volume_link();

            // Get incidence angle
            const vector3 normal = mask.normal(is.local);
            is.cos_incidence_angle = vector::dot(rd, normal);
        } else {
            is.status = false;
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
