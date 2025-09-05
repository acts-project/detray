/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/geometry/coordinates/cylindrical2D.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_cylinder_intersector.hpp"
#include "detray/tracks/helix.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <limits>
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t>
struct helix_intersector_impl;

/// @brief Intersection implementation for cylinder surfaces using helical
/// trajectories.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
/// @note Don't use for low p_t tracks!
template <algebra::concepts::aos algebra_t>
struct helix_intersector_impl<cylindrical2D<algebra_t>, algebra_t>
    : public ray_intersector_impl<cylindrical2D<algebra_t>, algebra_t,
                                  intersection::contains_pos> {

    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    template <typename surface_descr_t>
    using intersection_type =
        intersection2D<surface_descr_t, algebra_t, intersection::contains_pos>;

    template <typename other_algebra_t>
    using trajectory_type = detail::helix<other_algebra_t>;

    // Maximum number of solutions this intersector can produce
    static constexpr std::uint8_t n_solutions{2u};

    using result_type = darray<intersection_point_err<algebra_t>, n_solutions>;

    /// Operator function to find intersections between ray and planar mask
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tol negative cutoff for the path
    ///
    /// @return the intersection
    template <typename mask_t>
    DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
        const trajectory_type<algebra_t> &h, const transform3_type &trf,
        const mask_t &mask, const scalar_type = 0.f) const {

        result_type ret{};

        if (!run_rtsafe) {
            // Get the surface placement
            const auto &sm = trf.matrix();
            // Cylinder z axis
            const vector3_type sz = getter::vector<3>(sm, 0u, 2u);
            // Cylinder centre
            const point3_type sc = getter::vector<3>(sm, 0u, 3u);

            // Starting point on the helix for the Newton iteration
            // The mask is a cylinder -> it provides its radius as the first
            // value
            const scalar_type r{mask[cylinder2D::e_r]};

            // Try to guess the best starting positions for the iteration

            // Direction of the track at the helix origin
            const auto h_dir = h.dir(0.f);
            // Default starting path length for the Newton iteration (assumes
            // concentric cylinder)
            const scalar_type default_s{r * vector::perp(h_dir)};

            // Initial helix path length parameter
            darray<scalar_type, 2> paths{default_s, default_s};

            // try to guess good starting path by calculating the intersection
            // path of the helix tangential with the cylinder. This only has a
            // chance of working for tracks with reasonably high p_T !
            detail::ray<algebra_t> t{h.pos(), h.time(), h_dir, h.qop()};
            const auto qe = this->solve_intersection(t, mask, trf);

            // Obtain both possible solutions by looping over the (different)
            // starting positions
            auto n_runs{static_cast<unsigned int>(qe.solutions())};

            // Note: the default path length might be smaller than either
            // solution
            switch (qe.solutions()) {
                case 2:
                    paths[1] = qe.larger();
                    // If there are two solutions, reuse the case for a single
                    // solution to setup the intersection with the smaller path
                    // in ret[0]
                    [[fallthrough]];
                case 1: {
                    paths[0] = qe.smaller();
                    break;
                }
                    // Even if the ray is parallel to the cylinder, the helix
                    // might still hit it
                default: {
                    n_runs = 2u;
                    paths[0] = r;
                    paths[1] = -r;
                }
            }

            for (unsigned int i = 0u; i < n_runs; ++i) {

                scalar_type &s = paths[i];

                // Path length in the previous iteration step
                scalar_type s_prev{0.f};

                // f(s) = ((h.pos(s) - sc) x sz)^2 - r^2 == 0
                // Run the iteration on s
                std::size_t n_tries{0u};
                while (math::fabs(s - s_prev) > convergence_tolerance &&
                       n_tries < max_n_tries) {

                    // f'(s) = 2 * ( (h.pos(s) - sc) x sz) * (h.dir(s) x sz) )
                    const vector3_type crp = vector::cross(h.pos(s) - sc, sz);
                    const scalar_type denom{
                        2.f * vector::dot(crp, vector::cross(h.dir(s), sz))};

                    // No intersection can be found if dividing by zero
                    if (denom == 0.f) {
                        return {};
                    }

                    // x_n+1 = x_n - f(s) / f'(s)
                    s_prev = s;
                    s -= (vector::dot(crp, crp) - r * r) / denom;

                    ++n_tries;
                }
                // No intersection found within max number of trials
                if (n_tries == max_n_tries) {
                    return {};
                }

                ret[i] = {s, h.pos(s), s - s_prev};
            }

            return ret;
        } else {
            // Cylinder z axis
            const vector3_type sz = trf.z();
            // Cylinder centre
            const point3_type sc = trf.translation();

            // Starting point on the helix for the Newton iteration
            // The mask is a cylinder -> it provides its radius as the first
            // value
            const scalar_type r{mask[cylinder2D::e_r]};

            // Try to guess the best starting positions for the iteration

            // Direction of the track at the helix origin
            const auto h_dir = h.dir(0.5f * r);
            // Default starting path length for the Newton iteration (assumes
            // concentric cylinder)
            const scalar_type default_s{r * vector::perp(h_dir)};

            // Initial helix path length parameter
            darray<scalar_type, 2> paths{default_s, default_s};

            // try to guess good starting path by calculating the intersection
            // path of the helix tangential with the cylinder. This only has a
            // chance of working for tracks with reasonably high p_T !
            detail::ray<algebra_t> t{h.pos(), h.time(), h_dir, h.qop()};
            const auto qe = this->solve_intersection(t, mask, trf);

            // Obtain both possible solutions by looping over the (different)
            // starting positions
            auto n_runs{static_cast<unsigned int>(qe.solutions())};

            // Note: the default path length might be smaller than either
            // solution
            switch (qe.solutions()) {
                case 2:
                    paths[1] = qe.larger();
                    // If there are two solutions, reuse the case for a single
                    // solution to setup the intersection with the smaller path
                    // in ret[0]
                    [[fallthrough]];
                case 1: {
                    paths[0] = qe.smaller();
                    break;
                }
                default: {
                    n_runs = 2u;
                    paths[0] = r;
                    paths[1] = -r;
                }
            }

            /// Evaluate the function and its derivative at the point @param x
            auto cyl_inters_func = [&h, &r, &sz, &sc](const scalar_type x) {
                const vector3_type crp = vector::cross(h.pos(x) - sc, sz);

                // f(s) = ((h.pos(s) - sc) x sz)^2 - r^2 == 0
                const scalar_type f_s{(vector::dot(crp, crp) - r * r)};
                // f'(s) = 2 * ( (h.pos(s) - sc) x sz) * (h.dir(s) x sz) )
                const scalar_type df_s{
                    2.f * vector::dot(crp, vector::cross(h.dir(x), sz))};

                return std::make_tuple(f_s, df_s);
            };

            for (unsigned int i = 0u; i < n_runs; ++i) {

                const scalar_type &s_ini = paths[i];

                // Run the root finding algorithm
                const auto [s, ds] = newton_raphson_safe(cyl_inters_func, s_ini,
                                                         convergence_tolerance,
                                                         max_n_tries, max_path);

                ret[i] = {s, h.pos(s), ds};
            }

            return ret;
        }
    }

    /// Tolerance for convergence
    scalar_type convergence_tolerance{1.f * unit<scalar_type>::um};
    // Guard against inifinite loops
    std::size_t max_n_tries{1000u};
    // Early exit, if the intersection is too far away
    scalar_type max_path{5.f * unit<scalar_type>::m};
    // Complement the Newton algorithm with Bisection steps
    bool run_rtsafe{true};
};

template <concepts::algebra algebra_t>
struct helix_intersector_impl<concentric_cylindrical2D<algebra_t>, algebra_t>
    : public helix_intersector_impl<cylindrical2D<algebra_t>, algebra_t> {

    using base_type =
        helix_intersector_impl<cylindrical2D<algebra_t>, algebra_t>;

    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    template <typename surface_descr_t>
    using intersection_type =
        intersection2D<surface_descr_t, algebra_t, intersection::contains_pos>;

    template <typename other_algebra_t>
    using trajectory_type = detail::helix<other_algebra_t>;

    // Maximum number of solutions this intersector can produce
    static constexpr std::uint8_t n_solutions{1u};

    using result_type = intersection_point_err<algebra_t>;

    /// Operator function to find intersections between ray and planar mask
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tol negative cutoff for the path
    ///
    /// @return the intersection
    template <typename mask_t>
    DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
        const trajectory_type<algebra_t> &h, const transform3_type &trf,
        const mask_t &mask, const scalar_type = 0.f) const {

        // Array of two solutions
        typename base_type::result_type results =
            base_type::point_of_intersection(h, trf, mask, 0.f);

        // For portals, only take the solution along the helix direction,
        // not the one in the opposite direction
        return math::signbit(results[0].path) ? results[1] : results[0];
    }
};

}  // namespace detray
