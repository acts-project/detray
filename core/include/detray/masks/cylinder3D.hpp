/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical3.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_shape.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <string>
#include <tuple>

namespace detray {

/// @brief Geometrical shape of a full 3D cylinder.
///
/// It is defined by r and the two half lengths rel to the coordinate center.
class cylinder3D {
    public:
    /// The name for this shape
    inline static const std::string name = "cylinder3D";

    enum boundaries : std::size_t {
        e_r = 0,
        e_n_half_z = 1,
        e_p_half_z = 2,
        e_size = 3,
    };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = cylindrical3<algebra_t>;
    /// Local point type (3D)
    template <typename algebra_t>
    using loc_point_type = typename local_frame_type<algebra_t>::point3;

    /// Measurement frame
    template <typename algebra_t>
    using measurement_frame_type = local_frame_type<algebra_t>;
    /// Local measurement point (2D)
    template <typename algebra_t>
    using measurement_point_type = loc_point_type<algebra_t>;

    /// Underlying surface geometry: not a surface.
    template <typename algebra_t>
    using intersector_type = void;

    /// Behaviour of the three local axes (linear in r, circular in phi,
    /// linear in z)
    template <
        n_axis::shape e_s = n_axis::shape::e_open,
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular,
        template <typename, typename> class binning_loc2 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_r;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_phi;
        static constexpr n_axis::label axis_loc2 = n_axis::label::e_z;

        using types = std::tuple<n_axis::shape_t<e_s, axis_loc0>,
                                 n_axis::circular<axis_loc1>,
                                 n_axis::shape_t<e_s, axis_loc2>>;

        /// How to convert into the local axis system and back
        template <typename algebra_t>
        using coordinate_type = local_frame_type<algebra_t>;

        template <typename C, typename S>
        using binning = std::tuple<binning_loc0<C, S>, binning_loc1<C, S>,
                                   binning_loc2<C, S>>;
    };

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
        return (loc_p[0] <= bounds[e_r] + tol and
                bounds[e_n_half_z] - tol <= loc_p[2] and
                loc_p[2] <= bounds[e_p_half_z] + tol);
    }
};

}  // namespace detray
