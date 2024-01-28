/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/concentric_cylindrical2.hpp"
#include "detray/coordinates/cylindrical2.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/masks/cylinder2D.hpp"
#include "detray/navigation/detail/trajectories.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_cylinder_intersector.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <limits>
#include <type_traits>

namespace detray {

template <typename algebra_t, typename fame_t>
struct helix_intersector_impl;

/// @brief Intersection implementation for cylinder surfaces using helical
/// trajectories.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
/// @note Don't use for low p_t tracks!
template <typename algebra_t>
struct helix_intersector_impl<algebra_t, cylindrical2<algebra_t>>
    : public ray_intersector_impl<algebra_t, cylindrical2<algebra_t>> {

    using transform3_type = algebra_t;
    using scalar_type = typename transform3_type::scalar_type;
    using point2 = typename transform3_type::point2;
    using point3 = typename transform3_type::point3;
    using vector3 = typename transform3_type::vector3;

    template <typename surface_descr_t>
    using intersection_type = intersection2D<surface_descr_t, algebra_t>;
    using helix_type = detail::helix<transform3_type>;

    /// Operator function to find intersections between helix and cylinder mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_desc_t is the input transform type
    ///
    /// @param h is the input helix trajectory
    /// @param sf_desc is the surface descriptor
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    ///
    /// @return the intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline std::array<intersection_type<surface_descr_t>, 2>
    operator()(const helix_type &h, const surface_descr_t &sf_desc,
               const mask_t &mask, const transform3_type &trf,
               const scalar_type mask_tolerance = 0.f,
               const scalar_type = 0.f) const {

        std::array<intersection_type<surface_descr_t>, 2> ret;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{1000u};

        // Get the surface placement
        const auto &sm = trf.matrix();
        // Cylinder z axis
        const vector3 sz = getter::vector<3>(sm, 0u, 2u);
        // Cylinder centre
        const point3 sc = getter::vector<3>(sm, 0u, 3u);

        // Starting point on the helix for the Newton iteration
        // The mask is a cylinder -> it provides its radius as the first value
        const scalar_type r{mask[cylinder2D::e_r]};

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
        unsigned int n_runs =
            math::abs(paths[0] - paths[1]) < convergence_tolerance ? 1u : 2u;
        for (unsigned int i = 0u; i < n_runs; ++i) {

            scalar_type &s = paths[i];
            intersection_type<surface_descr_t> &is = ret[i];

            // Path length in the previous iteration step
            scalar_type s_prev{0.f};

            // f(s) = ((h.pos(s) - sc) x sz)^2 - r^2 == 0
            // Run the iteration on s
            std::size_t n_tries{0u};
            while (math::abs(s - s_prev) > convergence_tolerance and
                   n_tries < max_n_tries) {

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

            // Compute some additional information if the intersection is valid
            if (is.status == intersection::status::e_inside) {
                is.sf_desc = sf_desc;
                is.direction = math::signbit(s)
                                   ? intersection::direction::e_opposite
                                   : intersection::direction::e_along;
                is.volume_link = mask.volume_link();
            }
        }

        return ret;
    }

    /// Tolerance for convergence
    scalar_type convergence_tolerance{1e-3f};
};

template <typename algebra_t>
struct helix_intersector_impl<algebra_t, concentric_cylindrical2<algebra_t>>
    : public helix_intersector_impl<algebra_t, cylindrical2<algebra_t>> {};

}  // namespace detray
