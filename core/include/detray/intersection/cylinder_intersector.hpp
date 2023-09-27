/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical2.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// A functor to find intersections between a ray and a 2D cylinder mask
template <typename intersection_t>
struct cylinder_intersector {

    /// linear algebra types
    /// @{
    using transform3_type = typename intersection_t::transform3_type;
    using scalar_type = typename transform3_type::scalar_type;
    using point3 = typename transform3_type::point3;
    using point2 = typename transform3_type::point2;
    using vector3 = typename transform3_type::vector3;
    /// @}

    using intersection_type = intersection_t;
    using ray_type = detail::ray<transform3_type>;
    using helix_type = detail::helix<transform3_type>;

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
    template <typename mask_t, typename surface_t>
    DETRAY_HOST_DEVICE inline std::array<intersection_t, 2> operator()(
        const ray_type &ray, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        // One or both of these solutions might be invalid
        const auto qe = solve_intersection(ray, mask, trf);

        std::array<intersection_t, 2> ret;
        switch (qe.solutions()) {
            case 2:
                ret[1] = build_candidate(ray, mask, trf, qe.larger(),
                                         mask_tolerance);
                ret[1].sf_desc = sf;
                // If there are two solutions, reuse the case for a single
                // solution to setup the intersection with the smaller path
                // in ret[0]
                [[fallthrough]];
            case 1:
                ret[0] = build_candidate(ray, mask, trf, qe.smaller(),
                                         mask_tolerance);
                ret[0].sf_desc = sf;
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

    /// Operator function to find intersections between helix and cylinder mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_t is the input transform type
    ///
    /// @param h is the input helix trajectory
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    ///
    /// @return the intersection
    template <typename mask_t, typename surface_t>
    DETRAY_HOST_DEVICE inline std::array<intersection_t, 2> operator()(
        const helix_type &h, const surface_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        std::array<intersection_t, 2> ret;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{1000u};
        // Tolerance for convergence
        constexpr scalar_type tol{1e-4f};

        // Get the surface placement
        const auto &sm = trf.matrix();
        // Cylinder z axis
        const vector3 sz = getter::vector<3>(sm, 0u, 2u);
        // Cylinder centre
        const point3 sc = getter::vector<3>(sm, 0u, 3u);

        // Starting point on the helix for the Newton iteration
        // The mask is a cylinder -> it provides its radius as the first value
        const scalar_type r{mask[mask_t::boundaries::e_r]};

        // Try to guess the best starting positions for the iteration

        // Direction of the track at the helix origin
        const auto h_dir = h.dir(0.f);
        // Default starting path length for the Newton iteration (assumes
        // concentric cylinder)
        const scalar_type default_s{r * getter::perp(h_dir)};

        // Initial helix path length parameter
        std::array<scalar_type, 2> paths{default_s, default_s};

        // try to guess good starting path by calculating the intersection path
        // of the helix tangential with the cylinder. This only has a chance
        // of working for tracks with reasonably high p_T !
        detail::ray<transform3_type> t{h.pos(), h.time(), h_dir, h.qop()};
        const auto qe = this->solve_intersection(t, mask, trf);

        // Note: the default path length might be smaller than either solution
        switch (qe.solutions()) {
            case 2:
                paths[1] = qe.larger();
                // If there are two solutions, reuse the case for a single
                // solution to setup the intersection with the smaller path
                // in ret[0]
                [[fallthrough]];
            case 1:
                paths[0] = qe.smaller();
        };

        // Obtain both possible solutions by looping over the (different)
        // starting positions
        unsigned int n_runs = std::abs(paths[0] - paths[1]) < tol ? 1u : 2u;
        for (unsigned int i = 0u; i < n_runs; ++i) {

            scalar_type &s = paths[i];
            intersection_t &is = ret[i];

            // Path length in the previous iteration step
            scalar_type s_prev{0.f};

            // f(s) = ((h.pos(s) - sc) x sz)^2 - r^2 == 0
            // Run the iteration on s
            std::size_t n_tries{0u};
            while (std::abs(s - s_prev) > tol and n_tries < max_n_tries) {

                // f'(s) = 2 * ( (h.pos(s) - sc) x sz) * (h.dir(s) x sz) )
                const vector3 crp = vector::cross(h.pos(s) - sc, sz);
                const scalar_type denom{
                    2.f * vector::dot(crp, vector::cross(h.dir(s), sz))};
                // No intersection can be found if dividing by zero
                if (denom == 0.f) {
                    return ret;
                }
                // x_n+1 = x_n - f(s) / f'(s)
                s_prev = s;
                s -= (vector::dot(crp, crp) - r * r) / denom;

                ++n_tries;
            }
            // No intersection found within max number of trials
            if (n_tries == max_n_tries) {
                return ret;
            }

            // Build intersection struct from helix parameters
            is.path = s;
            const auto p3 = h.pos(s);
            is.local = mask.to_local_frame(trf, p3);
            is.status = mask.is_inside(is.local, mask_tolerance);

            // Perform the r-check for Newton solution even if it is not
            // required by the mask's shape
            const bool r_check =
                std::abs(r - is.local[2]) <
                mask_tolerance +
                    5.f * std::numeric_limits<scalar_type>::epsilon();
            if (not r_check) {
                is.status = intersection::status::e_outside;
            }

            // Compute some additional information if the intersection is valid
            if (is.status == intersection::status::e_inside) {
                is.sf_desc = sf;
                is.direction = std::signbit(s)
                                   ? intersection::direction::e_opposite
                                   : intersection::direction::e_along;
                is.volume_link = mask.volume_link();
            }
        }

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
              std::enable_if_t<std::is_same_v<typename mask_t::local_frame_type,
                                              cylindrical2<transform3_type>>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline void update(
        const ray_type &ray, intersection_t &sfi, const mask_t &mask,
        const transform3_type &trf,
        const scalar_type mask_tolerance = 0.f) const {

        // One or both of these solutions might be invalid
        const auto qe = solve_intersection(ray, mask, trf);

        switch (qe.solutions()) {
            case 1:
                sfi = build_candidate(ray, mask, trf, qe.smaller(),
                                      mask_tolerance);
                break;
            case 0:
                sfi.status = intersection::status::e_missed;
        };
    }

    protected:
    /// Calculates the distance to the (two) intersection points on the
    /// cylinder in global coordinates.
    ///
    /// @returns a quadratic equation object that contains the solution(s).
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline detail::quadratic_equation<scalar_type>
    solve_intersection(const ray_type &ray, const mask_t &mask,
                       const transform3_type &trf) const {
        const scalar_type r{mask[mask_t::shape::e_r]};
        const auto &m = trf.matrix();
        const vector3 sz = getter::vector<3>(m, 0u, 2u);
        const vector3 sc = getter::vector<3>(m, 0u, 3u);

        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();

        const auto pc_cross_sz = vector::cross(ro - sc, sz);
        const auto rd_cross_sz = vector::cross(rd, sz);
        const scalar_type a{vector::dot(rd_cross_sz, rd_cross_sz)};
        const scalar_type b{2.f * vector::dot(rd_cross_sz, pc_cross_sz)};
        const scalar_type c{vector::dot(pc_cross_sz, pc_cross_sz) - (r * r)};

        return detail::quadratic_equation<scalar_type>{a, b, c};
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

        // Construct the candidate only when needed
        if (path >= ray.overstep_tolerance()) {

            const point3 &ro = ray.pos();
            const vector3 &rd = ray.dir();

            is.path = path;
            const point3 p3 = ro + is.path * rd;

            is.local = mask.to_local_frame(trf, p3);
            is.status = mask.is_inside(is.local, mask_tolerance);

            // prepare some additional information in case the intersection
            // is valid
            if (is.status == intersection::status::e_inside) {
                is.direction = detail::signbit(is.path)
                                   ? intersection::direction::e_opposite
                                   : intersection::direction::e_along;
                is.volume_link = mask.volume_link();

                // Get incidence angle
                const scalar_type phi{is.local[0] / is.local[2]};
                const vector3 normal = {math_ns::cos(phi), math_ns::sin(phi),
                                        0.f};
                is.cos_incidence_angle = vector::dot(rd, normal);
            }
        } else {
            is.status = intersection::status::e_missed;
        }

        return is;
    }
};

}  // namespace detray