/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/coordinates/cartesian2D.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string>

namespace detray {

/// @brief Geometrical shape of a trapezoid2D.
///
/// It is defined by half lengths in local0 coordinate bounds[0] and bounds[1]
/// at -/+ half length in the local1 coordinate bounds[2]. bounds[3] contains
/// the precomputed value of 1 / (2 * bounds[2]), which avoids
/// excessive floating point divisions.
class trapezoid2D {
    public:
    /// The name for this shape
    inline static const std::string name = "trapezoid2D";

    enum boundaries : unsigned int {
        e_half_length_0 = 0u,
        e_half_length_1 = 1u,
        e_half_length_2 = 2u,
        e_divisor = 3u,  // 1 / (2 * bounds[e_half_length_2])
        e_size = 4u,
    };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = cartesian2D<algebra_t>;

    /// Dimension of the local coordinate system
    static constexpr std::size_t dim{2u};

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
    DETRAY_HOST_DEVICE inline auto check_boundaries(
        const bounds_t<scalar_t, kDIM> &bounds, const point_t &loc_p,
        const scalar_t tol = std::numeric_limits<scalar_t>::epsilon()) const {
        const scalar_t rel_y =
            (bounds[e_half_length_2] + loc_p[1]) * bounds[e_divisor];
        return (math::fabs(loc_p[0]) <= (bounds[e_half_length_0] +
                                         rel_y * (bounds[e_half_length_1] -
                                                  bounds[e_half_length_0]) +
                                         tol) &&
                math::fabs(loc_p[1]) <= (bounds[e_half_length_2] + tol));
    }

    /// @brief Measure of the shape: Area
    ///
    /// @param bounds the boundary values for this shape
    ///
    /// @returns the trapezoid area on the plane
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE constexpr scalar_t measure(
        const bounds_t<scalar_t, kDIM> &bounds) const {
        return area(bounds);
    }

    /// @brief The area of a the shape
    ///
    /// @param bounds the boundary values for this shape
    ///
    /// @returns the trapezoid area.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE constexpr scalar_t area(
        const bounds_t<scalar_t, kDIM> &bounds) const {
        return 2.f * (bounds[e_half_length_0] + bounds[e_half_length_1]) *
               bounds[e_half_length_2];
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
        const scalar_t x_bound{
            (bounds[e_half_length_0] > bounds[e_half_length_1]
                 ? bounds[e_half_length_0]
                 : bounds[e_half_length_1]) +
            env};
        const scalar_t y_bound{bounds[e_half_length_2] + env};
        return {-x_bound, -y_bound, -env, x_bound, y_bound, env};
    }

    /// @returns the shapes centroid in local cartesian coordinates
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE dpoint3D<algebra_t> centroid(
        const bounds_t<scalar_t, kDIM> &bounds) const {

        const scalar_t h_2{bounds[e_half_length_2]};
        const scalar_t a_2{bounds[e_half_length_1]};
        const scalar_t b_2{bounds[e_half_length_0]};

        const scalar_t y{2.f * h_2 * (2.f * a_2 + b_2) * 1.f /
                         (3.f * (a_2 + b_2))};

        return {0.f, y - h_2, 0.f};
    }

    /// Generate vertices in local cartesian frame
    ///
    /// @param bounds the boundary values for the trapezoid
    /// @param ls is the number of line segments
    ///
    /// @return a generated list of vertices
    template <typename point2_t, typename point3_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST dvector<point3_t> vertices(
        const bounds_t<scalar_t, kDIM> &bounds, dindex /*ignored*/) const {
        // left hand lower corner
        point3_t lh_lc{-bounds[e_half_length_0], -bounds[e_half_length_2], 0.f};
        // right hand lower corner
        point3_t rh_lc{bounds[e_half_length_0], -bounds[e_half_length_2], 0.f};
        // right hand upper corner
        point3_t rh_uc{bounds[e_half_length_1], bounds[e_half_length_2], 0.f};
        // left hand upper corner
        point3_t lh_uc{-bounds[e_half_length_1], bounds[e_half_length_2], 0.f};
        // Return the confining vertices
        return {lh_lc, rh_lc, rh_uc, lh_uc};
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

        constexpr auto tol{10.f * std::numeric_limits<scalar_t>::epsilon()};

        if (bounds[e_half_length_0] < tol or bounds[e_half_length_1] < tol) {
            os << "ERROR: Half length in x must be in the range (0, "
                  "numeric_max)"
               << std::endl;
            return false;
        }
        if (bounds[e_half_length_2] < tol) {
            os << "ERROR: Half length in y must be in the range (0, "
                  "numeric_max)"
               << std::endl;
            return false;
        }
        const auto div{1.f / (2.f * bounds[e_half_length_2])};
        if (math::fabs(bounds[e_divisor] - div) > tol) {
            os << "ERROR: Divisor incorrect. Should be: " << div << std::endl;
            return false;
        }

        return true;
    }
};

}  // namespace detray
