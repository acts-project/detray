/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/line2.hpp"
#include "detray/coordinates/coordinates.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/line_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_shape.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>

namespace detray {

/// @brief Mask for a line defined with a line length and its cross section
///
/// @tparam kSquareCrossSect determines whether the line has a cricular or
///         square cross section. This also changes the local coord. frame.
/// @tparam intersector_t defines how to intersect the underlying surface
///         geometry
///
/// The line can either have a circular or a square cross section. In the first
/// case bounds[0] refers to the radius, while in the second case it is the
/// half length of the square. The second boundary bounds[1] is the half length
/// in z.
template <bool kSquareCrossSect = false,
          template <typename> class intersector_t = line_intersector>
class line {
    public:
    /// The name for this shape
    inline static const std::string name = "line";

    /// Geometrical cross section of the line
    static constexpr bool square_cross_sect = kSquareCrossSect;

    enum boundaries : std::size_t {
        e_cross_section = 0,
        e_half_z = 1,
        e_size = 2
    };

    /// How to convert into the local system and back
    template <typename algebra_t>
    using local_frame_type = std::conditional_t<kSquareCrossSect, cartesian3<algebra_t>, line2<algebra_t>>;
    /// Measurement frame
    template <typename algebra_t>
    using measurement_frame_type = line2<algebra_t>;
    /// Local point type (2D)
    template <typename algebra_t>
    using loc_point_type = std::conditional_t<kSquareCrossSect, typename local_frame_type<algebra_t>::point3, typename local_frame_type<algebra_t>::point2>;
    /// Underlying surface geometry: planar
    template <typename algebra_t>
    using intersector_type = intersector_t<algebra_t>;

    /// Behaviour of the two local axes (linear in r/x, linear in z)
    template <
        n_axis::shape e_s = n_axis::shape::e_open,
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_r;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_z;

        using types = std::tuple<n_axis::shape_t<e_s, axis_loc0>,
                                 n_axis::shape_t<e_s, axis_loc1>>;

        /// How to convert into the local system and back
        template <typename algebra_t>
        using local_frame_type = std::conditional_t<kSquareCrossSect, cartesian3<algebra_t>, line2<algebra_t>>;

        template <typename C, typename S>
        using binning = std::tuple<binning_loc0<C, S>, binning_loc1<C, S>>;
    };

    /// @brief Check boundary values for a local point.
    ///
    /// @note the point is expected to be given in local coordinates by the
    /// caller. For the conversion from global cartesian coordinates, the
    /// nested @c shape struct can be used. In case of the line intersector
    /// this is the point of closest approach to the line.
    ///
    /// @tparam point_t the local point dimension can be 2 for a circular
    ///         cross section (local cylinder) or 3 for a square cross
    ///         section (local 3D cartesian).
    ///
    /// @param bounds the boundary values for this shape
    /// @param loc_p the point to be checked in the local coordinate system
    /// @param tol dynamic tolerance determined by caller
    ///
    /// @return true if the local point lies within the given boundaries.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM, typename point_t,
              typename std::enable_if_t<kDIM == 2, bool> = true>
    DETRAY_HOST_DEVICE inline bool check_boundaries(
        const bounds_t<scalar_t, kDIM> &bounds, const point_t &loc_p,
        const scalar_t tol = std::numeric_limits<scalar_t>::epsilon()) const {

        // For square cross section, we check if (1) x and y of the local cart.
        // point is less than the half cell size and (2) the distance to the
        // point of closest approach on thw line from the line center is less
        // than the half line length
        if constexpr (square_cross_sect) {
            return (std::abs(loc_p[0]) <= bounds[0] + tol &&
                    std::abs(loc_p[1]) <= bounds[0] + tol &&
                    std::abs(loc_p[2]) <= bounds[1] + tol);

            // For a circular cross section, we check if (1) the radial distance
            // is within the scope and (2) the distance to the point of closest
            // approach on the line from the line center is less than the line
            // half length
        } else {
            return (getter::perp(loc_p) <= bounds[0] + tol &&
                    std::abs(loc_p[2]) <= bounds[1] + tol);
        }
    }
};

}  // namespace detray