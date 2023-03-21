/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical3.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_bounds.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <string>

namespace detray {

/// @brief Geometrical shape of a full 3D cylinder.
///
/// It is defined by r and the two half lengths rel to the coordinate center.
class cylinder3D {
    public:
    /// The name for this shape
    inline static const std::string name = "cylinder3D";

    /// The measurement dimension (not allowed)
    inline static constexpr const unsigned int meas_dim{0u};

    enum boundaries : unsigned int {
        e_min_r = 0u,
        e_min_phi = 1u,
        e_min_z = 2u,
        e_max_r = 3u,
        e_max_phi = 4u,
        e_max_z = 5u,
        e_size = 6u,
    };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = cylindrical3<algebra_t>;
    /// Local point type (3D)
    template <typename algebra_t>
    using loc_point_type = typename local_frame_type<algebra_t>::point3;

    /// Measurement frame
    template <typename algebra_t>
    using measurement_frame_type = local_frame_type<algebra_t>;
    /// Local measurement point (2D)
    template <typename algebra_t>
    using measurement_point_type = loc_point_type<algebra_t>;

    /// Underlying surface geometry: not a surface.
    template <typename algebra_t>
    using intersector_type = void;

    /// Behaviour of the three local axes (linear in r, circular in phi,
    /// linear in z)
    template <
        n_axis::bounds e_s = n_axis::bounds::e_closed,
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular,
        template <typename, typename> class binning_loc2 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_r;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_phi;
        static constexpr n_axis::label axis_loc2 = n_axis::label::e_z;
        static constexpr std::size_t dim{3u};

        using types = dtuple<n_axis::bounds_t<e_s, axis_loc0>,
                             n_axis::circular<axis_loc1>,
                             n_axis::bounds_t<e_s, axis_loc2>>;

        /// How to convert into the local axis system and back
        template <typename algebra_t>
        using coordinate_type = local_frame_type<algebra_t>;

        template <typename C, typename S>
        using binning =
            dtuple<binning_loc0<C, S>, binning_loc1<C, S>, binning_loc2<C, S>>;
    };

    /// @brief Check boundary values for a local point.
    ///
    /// @note the point is expected to be given in local coordinates by the
    /// caller. For the conversion from global cartesian coordinates, the
    /// nested @c shape struct can be used. The point is assumed to be in
    /// the cylinder 3D frame.
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
        const bounds_t<scalar_t, kDIM> &bounds, const point_t &loc_p,
        const scalar_t tol = std::numeric_limits<scalar_t>::epsilon()) const {
        return (bounds[e_min_r] - tol <= loc_p[0] and
                bounds[e_min_phi] - tol <= loc_p[1] and
                bounds[e_min_z] - tol <= loc_p[2] and
                loc_p[0] <= bounds[e_max_r] + tol and
                loc_p[1] <= bounds[e_max_phi] + tol and
                loc_p[2] <= bounds[e_max_z] + tol);
    }

    /// @brief Lower and upper point for minimal axis aligned bounding box.
    ///
    /// Computes the min and max vertices in a local cartesian frame.
    ///
    /// @param bounds the boundary values for this shape
    /// @param env dynamic envelope around the shape
    ///
    /// @returns and array of coordinates that contains the lower point (first
    /// three values) and the upper point (latter three values).
    // @todo: Look at phi - range for a better fit
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE inline std::array<scalar_t, 6> local_min_bounds(
        const bounds_t<scalar_t, kDIM> &bounds,
        const scalar_t env = std::numeric_limits<scalar_t>::epsilon()) const {
        assert(env > 0.f);
        const scalar_t r_bound{bounds[e_max_r] + env};
        return {-r_bound, -r_bound, bounds[e_min_z] - env,
                r_bound,  r_bound,  bounds[e_max_z] + env};
    }
};

}  // namespace detray
