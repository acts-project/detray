/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cartesian2.hpp"
#include "detray/coordinates/cartesian3.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_bounds.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <ostream>
#include <string>

namespace detray {

/// @brief Underlying geometry for a single parameter bound mask
///
/// @tparam kCheckIndex is the index of the local point on which the mask is
///         applied
/// @tparam intersector_t defines how to intersect the underlying surface
///         geometry
template <unsigned int kCheckIndex = 0u,
          template <typename> class intersector_t = plane_intersector>
class single3D {
    public:
    /// The name for this shape
    inline static const std::string name = "single3D";

    enum boundaries : unsigned int {
        e_lower = 0u,
        e_upper = 1u,
        e_size = 2u,
    };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = cartesian2<algebra_t>;

    /// Underlying surface geometry: planar
    template <typename intersection_t>
    using intersector_type = intersector_t<intersection_t>;

    /// Behaviour of the two local axes (linear in single coordinate x, y or z)
    template <n_axis::bounds e_s = n_axis::bounds::e_closed,
              template <typename, typename> class binning_loc0 =
                  n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 =
            static_cast<n_axis::label>(kCheckIndex);
        static constexpr std::size_t dim{1u};

        /// How to convert into the local axis system and back
        template <typename algebra_t>
        using coordinate_type = local_frame_type<algebra_t>;

        using types = dtuple<n_axis::bounds_t<e_s, axis_loc0>>;

        template <typename C, typename S>
        using binning = dtuple<binning_loc0<C, S>>;
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
        const bounds_t<scalar_t, kDIM> &bounds, const point_t &loc_p,
        const scalar_t tol = std::numeric_limits<scalar_t>::epsilon()) const {
        return (bounds[e_lower] - tol <= loc_p[kCheckIndex] and
                loc_p[kCheckIndex] <= bounds[e_upper] + tol);
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
    DETRAY_HOST_DEVICE inline darray<scalar_t, 6> local_min_bounds(
        const bounds_t<scalar_t, kDIM> &bounds,
        const scalar_t env = std::numeric_limits<scalar_t>::epsilon()) const {
        assert(env > 0.f);
        darray<scalar_t, 6> o_bounds{-env, -env, -env, env, env, env};
        o_bounds[kCheckIndex] += bounds[e_lower];
        o_bounds[3u + kCheckIndex] += bounds[e_upper];
        return o_bounds;
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
        const bounds_t<scalar_t, kDIM> &bounds, std::ostream &os) const {

        if (bounds[e_upper] < bounds[e_lower]) {
            os << "ERROR: Upper bounds must be smaller than lower bounds ";
            return false;
        }

        return true;
    }
};

}  // namespace detray
