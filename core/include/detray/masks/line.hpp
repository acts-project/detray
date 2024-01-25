/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/cartesian3.hpp"
#include "detray/coordinates/line2.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/line_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_bounds.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <ostream>
#include <type_traits>

namespace detray {

/// @brief Geometrical shape of a line surface.
///
/// @tparam kSquareCrossSect determines whether the line has a cricular or
///         square cross section. This also changes the local coord. frame.
/// @tparam intersector_t defines how to intersect the underlying surface
///         geometry
/// @tparam kMeasDim defines the dimension of the measurement
/// @tparam kNormalOrder true if the index for measurement parameter follows
/// the local coordinate system
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

    enum boundaries : unsigned int {
        e_cross_section = 0u,
        e_half_z = 1u,
        e_size = 2u
    };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = line2<algebra_t>;

    /// Underlying surface geometry: line
    template <typename intersection_t>
    using intersector_type = intersector_t<intersection_t>;

    /// Behaviour of the two local axes (linear in r/x, linear in z)
    template <
        n_axis::bounds e_s = n_axis::bounds::e_closed,
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_r;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_z;
        static constexpr std::size_t dim{2u};

        using types = dtuple<n_axis::bounds_t<e_s, axis_loc0>,
                             n_axis::bounds_t<e_s, axis_loc1>>;

        /// How to convert into the local axis system and back
        template <typename algebra_t>
        using coordinate_type = line2<algebra_t>;

        template <typename C, typename S>
        using binning = dtuple<binning_loc0<C, S>, binning_loc1<C, S>>;
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
              typename std::enable_if_t<kDIM == 2u, bool> = true>
    DETRAY_HOST_DEVICE inline bool check_boundaries(
        const bounds_t<scalar_t, kDIM> &bounds, const point_t &loc_p,
        const scalar_t tol = std::numeric_limits<scalar_t>::epsilon()) const {

        // For a square cross section (e.g. a cell of drift chamber), we check
        // if (1) x and y of the local cart. point is less than the half cell
        // size and (2) the distance to the point of closest approach on thw
        // line from the line center is less than the half line length
        if constexpr (square_cross_sect) {
            return (math::abs(loc_p[0] * math::cos(loc_p[2])) <=
                        bounds[e_cross_section] + tol &&
                    math::abs(loc_p[0] * math::sin(loc_p[2])) <=
                        bounds[e_cross_section] + tol &&
                    math::abs(loc_p[1]) <= bounds[e_half_z] + tol);

        }
        // For a circular cross section (e.g. straw tube), we check if (1) the
        // radial distance is within the scope and (2) the distance to the point
        // of closest approach on the line from the line center is less than the
        // line half length
        else {
            return (loc_p[0] <= bounds[e_cross_section] + tol &&
                    math::abs(loc_p[1]) <= bounds[e_half_z] + tol);
        }
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
        const scalar_t xy_bound{bounds[e_cross_section] + env};
        const scalar_t z_bound{bounds[e_half_z] + env};
        return {-xy_bound, -xy_bound, -z_bound, xy_bound, xy_bound, z_bound};
    }

    /// @returns the shapes centroid in local cartesian coordinates
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE typename algebra_t::point3 centroid(
        const bounds_t<scalar_t, kDIM> &) const {

        return {0.f, 0.f, 0.f};
    }

    /// Generate vertices in local cartesian frame
    ///
    /// @param bounds the boundary values for the line
    /// @param n_seg is the number of line segments
    ///
    /// @return a generated list of vertices
    template <typename point2_t, typename point3_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST dvector<point3_t> vertices(
        const bounds_t<scalar_t, kDIM> &bounds, dindex /*ignored*/) const {

        point3_t lc = {0.f, 0.f, -bounds[e_half_z]};
        point3_t rc = {0.f, 0.f, bounds[e_half_z]};

        return {lc, rc};
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

        if (bounds[e_cross_section] < tol) {
            os << "ERROR: Radius/sides must be in the range (0, numeric_max)"
               << std::endl;
            return false;
        }
        if (bounds[e_half_z] < tol) {
            os << "ERROR: Half length z must be in the range (0, numeric_max)"
               << std::endl;
            return false;
        }

        return true;
    }
};

}  // namespace detray
