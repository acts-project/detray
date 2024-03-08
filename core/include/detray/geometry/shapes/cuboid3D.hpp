/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cartesian3.hpp"
#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string>

namespace detray {

/// @brief Geometrical shape of a full 3D cuboid.
///
/// It is defined by 3 min-max length pairs and checks whether a point is
/// somewhere inside the cuboid. This type is mainly used for aabb description.
class cuboid3D {
    public:
    /// The name for this shape
    inline static const std::string name = "cuboid3D";

    enum boundaries : unsigned int {
        e_min_x = 0u,
        e_min_y = 1u,
        e_min_z = 2u,
        e_max_x = 3u,
        e_max_y = 4u,
        e_max_z = 5u,
        e_size = 6u,
    };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = cartesian3<algebra_t>;

    /// Dimension of the local coordinate system
    static constexpr std::size_t dim{3u};

    /// @brief Check boundary values for a local point.
    ///
    /// @note the point is expected to be given in local coordinates by the
    /// caller. For the conversion from global cartesian coordinates, the
    /// nested @c shape struct can be used. The point is assumed to be in
    /// the a cartesian 3D frame.
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
        return (bounds[e_min_x] - tol <= loc_p[0] and
                bounds[e_min_y] - tol <= loc_p[1] and
                bounds[e_min_x] - tol <= loc_p[2] and
                loc_p[0] <= bounds[e_max_x] + tol and
                loc_p[1] <= bounds[e_max_y] + tol and
                loc_p[2] <= bounds[e_max_z] + tol);
    }

    /// @brief Measure of the shape: Volume
    ///
    /// @param bounds the boundary values for this shape
    ///
    /// @returns the cuboid volume as part of global space.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE constexpr scalar_t measure(
        const bounds_t<scalar_t, kDIM> &bounds) const {
        return volume(bounds);
    }

    /// @brief The volume of a the shape
    ///
    /// @param bounds the boundary values for this shape
    ///
    /// @returns the cuboid volume.
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE constexpr scalar_t volume(
        const bounds_t<scalar_t, kDIM> &bounds) const {
        return (bounds[e_max_x] - bounds[e_min_x]) *
               (bounds[e_max_y] - bounds[e_min_y]) *
               (bounds[e_max_z] - bounds[e_min_z]);
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
        bounds_t<scalar_t, kDIM> o_bounds{bounds};
        for (unsigned int i{0}; i < 3u; ++i) {
            o_bounds[i] -= env;
            o_bounds[i + 3u] += env;
        }
        return o_bounds;
    }

    /// @returns the shapes centroid in local cartesian coordinates
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE typename algebra_t::point3 centroid(
        const bounds_t<scalar_t, kDIM> &bounds) const {

        return 0.5f *
               typename algebra_t::point3{bounds[e_min_x] + bounds[e_max_x],
                                          bounds[e_min_y] + bounds[e_max_y],
                                          bounds[e_min_z] + bounds[e_max_z]};
    }

    /// Generate vertices in local cartesian frame
    ///
    /// @param bounds the boundary values for the cuboid
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
            "Vertex generation for cuboids is not implemented");
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

        if (bounds[e_min_x] >= bounds[e_max_x] or
            math::abs(bounds[e_min_x] - bounds[e_max_x]) < tol) {
            os << "ERROR: Min x must be smaller than max x.";
            return false;
        }
        if (bounds[e_min_y] >= bounds[e_max_y] or
            math::abs(bounds[e_min_y] - bounds[e_max_y]) < tol) {
            os << "ERROR: Min y must be smaller than max y.";
            return false;
        }
        if (bounds[e_min_z] >= bounds[e_max_z] or
            math::abs(bounds[e_min_z] - bounds[e_max_z]) < tol) {
            os << "ERROR: Min z must be smaller than max z.";
            return false;
        }

        return true;
    }
};

}  // namespace detray
