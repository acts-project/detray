/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical3.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/cylinder_intersector.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <ostream>
#include <string>

namespace detray {

/// @brief Geometrical shape of a full 3D cylinder.
///
/// It is defined by r and the two half lengths rel to the coordinate center.
class cylinder3D {
    public:
    /// The name for this shape
    inline static const std::string name = "cylinder3D";

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

    /// Underlying surface geometry: not a surface.
    template <typename intersection_t>
    using intersector_type = void;

    /// Dimension of the local coordinate system
    static constexpr std::size_t dim{3u};

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
    DETRAY_HOST_DEVICE inline darray<scalar_t, 6> local_min_bounds(
        const bounds_t<scalar_t, kDIM> &bounds,
        const scalar_t env = std::numeric_limits<scalar_t>::epsilon()) const {
        assert(env > 0.f);
        const scalar_t r_bound{bounds[e_max_r] + env};
        return {-r_bound, -r_bound, bounds[e_min_z] - env,
                r_bound,  r_bound,  bounds[e_max_z] + env};
    }

    /// @returns the shapes centroid in local cartesian coordinates
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE typename algebra_t::point3 centroid(
        const bounds_t<scalar_t, kDIM> &bounds) const {

        return 0.5f * typename algebra_t::point3{
                          0.f, (bounds[e_min_phi] + bounds[e_max_phi]),
                          (bounds[e_min_z] + bounds[e_max_z])};
    }

    /// Generate vertices in local cartesian frame
    ///
    /// @param bounds the boundary values for the cylinder
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
            "Vertex generation for 3D cylinders is not implemented");
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

        if (bounds[e_min_r] < tol) {
            os << "ERROR: Radii must be in the range (0, numeric_max)"
               << std::endl;
            return false;
        }
        if (bounds[e_min_r] >= bounds[e_max_r] or
            math::abs(bounds[e_min_r] - bounds[e_max_r]) < tol) {
            os << "ERROR: Min Radius must be smaller than max Radius.";
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
