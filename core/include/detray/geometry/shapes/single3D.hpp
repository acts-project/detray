/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cartesian2.hpp"
#include "detray/coordinates/cartesian3.hpp"
#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string>

namespace detray {

/// @brief Underlying geometry for a single parameter bound mask
///
/// @tparam kCheckIndex is the index of the local point on which the mask is
///         applied
template <unsigned int kCheckIndex = 0u>
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

    /// Dimension of the local coordinate system
    static constexpr std::size_t dim{1u};

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

    /// @returns the shapes centroid in local cartesian coordinates
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE auto centroid(
        const bounds_t<scalar_t, kDIM> &bounds) const {

        using point3_t = typename algebra_t::point3;

        point3_t centr{0.f, 0.f, 0.f};
        centr[kCheckIndex] = 0.5f * (bounds[e_lower] + bounds[e_upper]);

        return centr;
    }

    /// Generate vertices in local cartesian frame
    ///
    /// @param bounds the boundary values for the single value
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
            "Vertex generation for single value shapes is not implemented");
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

        if (bounds[e_upper] < bounds[e_lower]) {
            os << "ERROR: Upper bounds must be smaller than lower bounds ";
            return false;
        }

        return true;
    }
};

}  // namespace detray
