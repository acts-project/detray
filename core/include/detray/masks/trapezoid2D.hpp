/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cartesian2.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_bounds.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <string>
#include <tuple>

namespace detray {

/// @brief Geometrical shape of a trapezoid2D.
///
/// @tparam intersector_t defines how to intersect the underlying surface
///         geometry
/// @tparam kMeasDim defines the dimension of the measurement
///
/// It is defined by half lengths in local0 coordinate bounds[0] and bounds[1]
/// at -/+ half length in the local1 coordinate bounds[2]. bounds[3] contains
/// the precomputed value of 1 / (2 * bounds[2]), which avoids
/// excessive floating point divisions.
template <template <typename> class intersector_t = plane_intersector,
          std::size_t kMeasDim = 2>
class trapezoid2D {
    public:
    /// The name for this shape
    inline static const std::string name = "trapezoid2D";

    /// The measurement dimension
    inline static constexpr const std::size_t meas_dim = kMeasDim;

    enum boundaries : std::size_t {
        e_half_length_0 = 0,
        e_half_length_1 = 1,
        e_half_length_2 = 2,
        e_divisor = 3,  // 1 / (2 * bounds[e_half_length_2])
        e_size = 4,
    };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = cartesian2<algebra_t>;
    /// Local point type (2D)
    template <typename algebra_t>
    using loc_point_type = typename local_frame_type<algebra_t>::point2;

    /// Measurement frame
    template <typename algebra_t>
    using measurement_frame_type = local_frame_type<algebra_t>;
    /// Local measurement point (2D)
    template <typename algebra_t>
    using measurement_point_type = loc_point_type<algebra_t>;

    /// Underlying surface geometry: planar
    template <typename algebra_t>
    using intersector_type = intersector_t<algebra_t>;

    /// Behaviour of the two local axes (linear in x, y)
    template <
        n_axis::bounds e_s = n_axis::bounds::e_closed,
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_x;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_y;
        static constexpr std::size_t dim{2UL};

        /// How to convert into the local axis system and back
        template <typename algebra_t>
        using coordinate_type = local_frame_type<algebra_t>;

        using types = std::tuple<n_axis::bounds_t<e_s, axis_loc0>,
                                 n_axis::bounds_t<e_s, axis_loc1>>;

        template <typename C, typename S>
        using binning = std::tuple<binning_loc0<C, S>, binning_loc1<C, S>>;
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
        scalar_t rel_y =
            (bounds[e_half_length_2] + loc_p[1]) * bounds[e_divisor];
        return (std::abs(loc_p[0]) <= bounds[e_half_length_0] +
                                          rel_y * (bounds[e_half_length_1] -
                                                   bounds[e_half_length_0]) +
                                          tol and
                std::abs(loc_p[1]) <= bounds[e_half_length_2] + tol);
    }

    template <typename param_t>
    DETRAY_HOST_DEVICE inline typename param_t::point2 to_measurement(
        param_t& param, const typename param_t::point2& offset = {0, 0}) const {
        return param.local() + offset;
    }
};

}  // namespace detray
