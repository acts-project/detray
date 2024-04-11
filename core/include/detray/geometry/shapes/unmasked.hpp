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
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/coordinates/cartesian2D.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string>

namespace detray {

/// @brief Generic shape without boundaries.
template <std::size_t DIM = 2>
class unmasked {
    public:
    /// The name for this shape
    inline static const std::string name = "unmasked";

    enum boundaries : unsigned int { e_size = 1u };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = cartesian2D<algebra_t>;

    /// Dimension of the local coordinate system
    static constexpr std::size_t dim{DIM};

    /// @brief Check boundary values for a local point.
    ///
    /// @tparam bounds_t any type of boundary values
    ///
    /// @note the parameters are ignored
    ///
    /// @return always true
    template <typename bounds_t, typename point_t, typename scalar_t>
    DETRAY_HOST_DEVICE inline constexpr bool check_boundaries(
        const bounds_t& /*bounds*/, const point_t& /*loc_p*/,
        const scalar_t /*tol*/) const {
        return true;
    }

    /// @brief Measure of the shape: Inf
    ///
    /// @param bounds the boundary values for this shape
    ///
    /// @returns Inf.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM>
    DETRAY_HOST_DEVICE constexpr scalar_t measure(
        const bounds_t<scalar_t, kDIM>& bounds) const {
        if constexpr (dim == 2) {
            return area(bounds);
        } else {
            return volume(bounds);
        }
    }

    /// @brief The area of a the shape
    ///
    /// @param bounds the boundary values for this shape
    ///
    /// @returns the stereo annulus area.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM, std::size_t D = dim,
              std::enable_if_t<D == 2, bool> = true>
    DETRAY_HOST_DEVICE constexpr scalar_t area(
        const bounds_t<scalar_t, kDIM>&) const {
        return std::numeric_limits<scalar_t>::max();
    }

    /// @brief The volume of a the shape
    ///
    /// @param bounds the boundary values for this shape
    ///
    /// @returns Inf.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM, std::size_t D = dim,
              std::enable_if_t<D == 3, bool> = true>
    DETRAY_HOST_DEVICE constexpr scalar_t volume(
        const bounds_t<scalar_t, kDIM>&) const {
        return std::numeric_limits<scalar_t>::max();
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
    DETRAY_HOST_DEVICE constexpr darray<scalar_t, 6> local_min_bounds(
        const bounds_t<scalar_t, kDIM>& /*bounds*/,
        const scalar_t /*env*/ =
            std::numeric_limits<scalar_t>::epsilon()) const {
        constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
        return {-inv, -inv, -inv, inv, inv, inv};
    }

    /// @returns the shapes centroid in global cartesian coordinates
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE dpoint3D<algebra_t> centroid(
        const bounds_t<scalar_t, kDIM>&) const {
        return {0.f, 0.f, 0.f};
    }

    /// Generate vertices in local cartesian frame
    ///
    /// @param bounds the boundary values
    /// @param n_seg is the number of line segments
    ///
    /// @return a generated list of vertices
    template <typename point2_t, typename point3_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST dvector<point3_t> vertices(
        const bounds_t<scalar_t, kDIM>& bounds, dindex) const {
        return local_min_bounds(bounds);
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
        const bounds_t<scalar_t, kDIM>& /*bounds*/,
        std::ostream& /*os*/) const {
        return true;
    }
};

}  // namespace detray
