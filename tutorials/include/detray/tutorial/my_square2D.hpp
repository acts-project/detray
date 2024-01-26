/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cartesian3.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/plane_intersector.hpp"

// System include(s)
#include <limits>
#include <string>

namespace detray::tutorial {

/// @brief Geometrical shape of a square.
///
/// @tparam intersector_t defines how to intersect the underlying surface
///         geometry
/// @tparam kMeasDim defines the dimension of the measurement
/// @tparam kNormalOrder true if the index for measurement parameter follows
/// the local coordinate system
///
/// It is defined by the half length of the square (bounds[0]),
/// and can be checked with a tolerance in t.
template <template <typename> class intersector_t = plane_intersector,
          unsigned int kMeasDim = 2u>
class square2D {
    public:
    /// The name for this shape
    inline static const std::string name = "square2D";

    /// The measurement dimension
    inline static constexpr const unsigned int meas_dim{kMeasDim};

    /// Normal ordering
    inline static constexpr const bool normal_order{true};

    // Measurement dimension check
    static_assert(meas_dim == 1u || meas_dim == 2u,
                  "Only 1D or 2D measurement is allowed");

    enum boundaries : unsigned int {
        e_half_length = 0,  // < boundary value: the half length of the square
        e_size = 1u,        // < Number of boundary values for this shape
    };

    /// Local coordinate frame for boundary checks: cartesian
    template <typename algebra_t>
    using local_frame_type = cartesian3<algebra_t>;

    /// Underlying surface geometry: planar
    template <typename intersection_t>
    using intersector_type = intersector_t<intersection_t>;

    /// Behaviour of the two local axes (linear in x, linear in y)
    /// Needed to map a grid onto the square (material maps)
    template <
        n_axis::bounds e_s = n_axis::bounds::e_closed,
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_x;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_y;
        static constexpr std::size_t dim{2u};

        /// How to convert into the local axis system and back
        template <typename algebra_t>
        using coordinate_type = local_frame_type<algebra_t>;

        using types = dtuple<n_axis::bounds_t<e_s, axis_loc0>,
                             n_axis::bounds_t<e_s, axis_loc1>>;

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
        return (math::abs(loc_p[0]) <= bounds[e_half_length] + tol and
                math::abs(loc_p[1]) <= bounds[e_half_length] + tol);
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
        const scalar_t bound{bounds[e_half_length] + env};
        return {-bound, -bound, -env, bound, bound, env};
    }
};

}  // namespace detray::tutorial
