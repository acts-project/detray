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
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string>

namespace detray {

/// @brief Geometrical shape of a closed ring.
///
/// It is defined by the two radii bounds[0] and bounds[1],
/// and can be checked with a tolerance in t[0] and t[1].
class ring2D {
    public:
    /// The name for this shape
    inline static const std::string name = "ring2D";

    enum boundaries : unsigned int {
        e_inner_r = 0u,
        e_outer_r = 1u,
        e_size = 2u,
    };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = polar2D<algebra_t>;

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
    DETRAY_HOST_DEVICE inline bool check_boundaries(
        const bounds_t<scalar_t, kDIM> &bounds, const point_t &loc_p,
        const scalar_t tol = std::numeric_limits<scalar_t>::epsilon()) const {

        return (loc_p[0] + tol >= bounds[e_inner_r] and
                loc_p[0] <= bounds[e_outer_r] + tol);
    }

    /// @brief Measure of the shape: Area
    ///
    /// @param bounds the boundary values for this shape
    ///
    /// @returns the ring area on the plane
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
    /// @returns the ring area.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE constexpr scalar_t area(
        const bounds_t<scalar_t, kDIM> &bounds) const {
        return (bounds[e_outer_r] * bounds[e_outer_r] -
                bounds[e_inner_r] * bounds[e_inner_r]) *
               constant<scalar>::pi;
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
        const scalar_t r_bound{env + bounds[e_outer_r]};
        return {-r_bound, -r_bound, -env, r_bound, r_bound, env};
    }

    /// @returns the shapes centroid in local cartesian coordinates
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE dpoint3D<algebra_t> centroid(
        const bounds_t<scalar_t, kDIM> &) const {

        return {0.f, 0.f, 0.f};
    }

    /// Generate vertices in local cartesian frame
    ///
    /// @param bounds the boundary values for the stereo annulus
    /// @param n_seg is the number of line segments
    ///
    /// @return a generated list of vertices
    template <typename point2_t, typename point3_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST dvector<point3_t> vertices(const bounds_t<scalar_t, kDIM> &,
                                           dindex) const {
        throw std::runtime_error(
            "Vertex generation for rings is not implemented");
        return {};
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

        if (math::signbit(bounds[e_inner_r]) or bounds[e_outer_r] < tol) {
            os << "ERROR: Radius must be in the range [0, numeric_max)"
               << std::endl;
            return false;
        }
        if (bounds[e_inner_r] >= bounds[e_outer_r] or
            math::fabs(bounds[e_inner_r] - bounds[e_outer_r]) < tol) {
            os << "ERROR: Inner radius must be smaller outer radius.";
            return false;
        }

        return true;
    }
};

}  // namespace detray
