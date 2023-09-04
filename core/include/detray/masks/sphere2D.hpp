/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/spherical3D.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/sphere_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_bounds.hpp"

// System include(s)
#include <limits>
#include <string>

namespace detray {

/// @brief Geometrical shape of a 2 dimensional sphere.
///
/// @tparam intersector_t defines how to intersect the underlying surface
///         geometry
template <template <typename> class intersector_t = sphere_intersector>
class sphere2D {
    public:
    /// The name for this shape
    inline static const std::string name = "sphere2D";

    /// Names for the mask boundary values
    enum boundaries : unsigned int {
        e_r = 0u,
        e_size = 2u,
    };

    /// Local coordinate frame ( spherical )
    template <typename algebra_t>
    using local_frame_type = spherical3D<algebra_t>;

    /// Underlying surface geometry: spherical
    template <typename intersection_t>
    using intersector_type = intersector_t<intersection_t>;

    /// Behaviour of the two local axes (circular in phi, circular in theta)
    template <
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_phi;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_theta;
        static constexpr std::size_t dim{2u};

        using types =
            dtuple<n_axis::circular<axis_loc0>, n_axis::circular<axis_loc1>>;

        /// Local coordinate frame (both for disc and focal system ?)
        template <typename algebra_t>
        using coordinate_type = local_frame_type<algebra_t>;

        template <typename C, typename S>
        using binning = dtuple<binning_loc0<C, S>, binning_loc1<C, S>>;
    };

    /// @brief Check boundary values for a local point.
    ///
    /// @note the point is expected to be given in local coordinates by the
    /// caller. For the conversion from global cartesian coordinates, the
    /// nested @c shape struct can be used.
    ///
    /// @param bounds the boundary values for this shape
    /// @param loc_p the point to be checked in the local coordinate system
    /// @param tol dynamic tolerance determined by caller
    ///
    /// @return true if the local point lies within the given boundaries.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM, typename point_t,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE inline bool check_boundaries(
        const bounds_t<scalar_t, kDIM>& bounds, const point_t& loc_p,
        const scalar_t tol = std::numeric_limits<scalar_t>::epsilon()) const {
        return (math_ns::abs(loc_p[0] - bounds[e_r]) > tol);
    }

    /// @brief Lower and upper point for minimal axis aligned bounding box.
    ///
    /// Computes the min and max vertices in a local cartesian frame.
    ///
    /// @param bounds the boundary values for this shape
    /// @param env dynamic envelope around the shape
    ///
    /// @returns and array of coordinates that contains the lower point (first
    /// three values) and the upper point (latter three values) .
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE inline std::array<scalar_t, 6> local_min_bounds(
        const bounds_t<scalar_t, kDIM>& bounds,
        const scalar_t env = std::numeric_limits<scalar_t>::epsilon()) const {
        assert(env > 0.f);
        const scalar_t r_bound = env + bounds[e_r];
        return {-r_bound, -r_bound, -env, r_bound, r_bound, env};
    }

    /// @brief Check consistency of boundary values.
    ///
    /// @param bounds the boundary values for this shape
    /// @param os output stream for error messages
    ///
    /// @return true if the bounds are consistent.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST constexpr bool check_consistency(
        const bounds_t<scalar_t, kDIM>& bounds, std::ostream& os) const {

        constexpr auto tol{10.f * std::numeric_limits<scalar_t>::epsilon()};

        if (bounds[e_r] < tol) {
            os << "ERROR: Radius must be in the range [0, numeric_max)"
               << std::endl;
            return false;
        }

        return true;
    }
};

}  // namespace detray
