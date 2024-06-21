/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/navigation/detail/helix.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/utils/root_finding.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, typename algebra_t>
struct helix_intersector_impl;

/// @brief Intersection implementation for helical trajectories with planar
/// surfaces.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
template <typename algebra_t>
struct helix_intersector_impl<cartesian2D<algebra_t>, algebra_t> {

    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    template <typename surface_descr_t>
    using intersection_type = intersection2D<surface_descr_t, algebra_t>;
    using helix_type = detail::helix<algebra_t>;

    /// Operator function to find intersections between helix and planar mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_desc_t is the input surface descriptor type
    ///
    /// @param h is the input helix trajectory
    /// @param sf_desc is the surface descriptor
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tolerance is the tolerance for track overstepping
    ///
    /// @return the intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const helix_type &h, const surface_descr_t &sf_desc, const mask_t &mask,
        const transform3_type &trf,
        const std::array<scalar_type, 2u> mask_tolerance =
            {detail::invalid_value<scalar_type>(),
             detail::invalid_value<scalar_type>()},
        const scalar_type = 0.f, const scalar_type = 0.f) const {

        assert((mask_tolerance[0] == mask_tolerance[1]) &&
               "Helix intersectors use only one mask tolerance value");

        intersection_type<surface_descr_t> sfi;

        // Surface normal
        const vector3_type sn = trf.z();
        // Surface translation
        const point3_type st = trf.translation();

        // Starting point on the helix for the Newton iteration
        const vector3_type dist{trf.point_to_global(mask.centroid()) -
                                h.pos(0.f)};
        scalar_type denom{vector::dot(sn, h.dir(0.5f * getter::norm(dist)))};
        scalar_type s_ini;
        if (denom == 0.f) {
            s_ini = getter::norm(dist);
        } else {
            s_ini = vector::dot(sn, dist) / denom;
        }

        /// Evaluate the function and its derivative at the point @param x
        auto plane_inters_func = [&h, &st, &sn](const scalar_type x) {
            // f(s) = sn * (h.pos(s) - st) == 0
            const scalar_type f_s{vector::dot(sn, (h.pos(x) - st))};
            // f'(s) = sn * h.dir(s)
            const scalar_type df_s{vector::dot(sn, h.dir(x))};

            return std::make_tuple(f_s, df_s);
        };

        // Run the root finding algorithm
        const auto [s, ds] =
            newton_raphson_safe(plane_inters_func, s_ini, convergence_tolerance,
                                max_n_tries, max_path);

        // Build intersection struct from the root
        build_intersection(h, sfi, s, ds, sf_desc, mask, trf, mask_tolerance);

        return sfi;
    }

    /// Interface to use fixed mask tolerance
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const helix_type &h, const surface_descr_t &sf_desc, const mask_t &mask,
        const transform3_type &trf, const scalar_type mask_tolerance,
        const scalar_type = 0.f, const scalar_type = 0.f) const {
        return this->operator()(h, sf_desc, mask, trf,
                                {mask_tolerance, mask_tolerance}, 0.f);
    }

    /// Tolerance for convergence
    scalar_type convergence_tolerance{1.f * unit<scalar_type>::um};
    // Guard against inifinite loops
    std::size_t max_n_tries{1000u};
    // Early exit, if the intersection is too far away
    scalar_type max_path{5.f * unit<scalar_type>::m};
};

template <typename algebra_t>
struct helix_intersector_impl<polar2D<algebra_t>, algebra_t>
    : public helix_intersector_impl<cartesian2D<algebra_t>, algebra_t> {};

}  // namespace detray
