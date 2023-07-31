/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/cylindrical2.hpp"
#include "detray/coordinates/cylindrical3.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_bounds.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <ostream>
#include <string>

namespace detray {

/// @brief Geometrical shape of a 2D cylinder.
///
/// @tparam kRadialCheck is a boolean to steer whether the radius compatibility
///         needs to be checked (changes local coordinate system def.)
/// @tparam intersector_t defines how to intersect the underlying surface
///         geometry
/// @tparam kMeasDim defines the dimension of the measurement
/// @tparam kNormalOrder true if the index for measurement parameter follows
/// the local coordinate system
///
/// It is defined by r and the two half lengths rel to the coordinate center.
template <bool kRadialCheck = false,
          template <typename> class intersector_t = cylinder_intersector,
          unsigned int kMeasDim = 2u, bool kNormalOrder = true>
class cylinder2D {
    public:
    /// The name for this shape
    inline static const std::string name = "cylinder2D";

    /// Check the radial position in boundary check
    static constexpr bool check_radius = kRadialCheck;

    /// The measurement dimension
    inline static constexpr const unsigned int meas_dim{kMeasDim};

    /// Normal ordering
    inline static constexpr const bool normal_order{kNormalOrder};

    // Measurement dimension check
    static_assert(meas_dim == 1u || meas_dim == 2u,
                  "Only 1D or 2D measurement is allowed");

    enum boundaries : unsigned int {
        e_r = 0u,
        e_n_half_z = 1u,
        e_p_half_z = 2u,
        e_size = 3u,
    };

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type = cylindrical2<algebra_t>;

    /// Underlying surface geometry: cylindrical
    template <typename intersection_t>
    using intersector_type = intersector_t<intersection_t>;

    /// Behaviour of the two local axes (circular in r_phi, linear in z)
    template <
        n_axis::bounds e_s = n_axis::bounds::e_closed,
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_rphi;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_cyl_z;
        static constexpr std::size_t dim{2u};

        using types = dtuple<n_axis::circular<axis_loc0>,
                             n_axis::bounds_t<e_s, axis_loc1>>;

        /// How to convert into the local axis system and back
        template <typename algebra_t>
        using coordinate_type = cylindrical2<algebra_t>;

        template <typename C, typename S>
        using binning = dtuple<binning_loc0<C, S>, binning_loc1<C, S>>;
    };

    /// @brief Check boundary values for a local point.
    ///
    /// @note the point is expected to be given in local coordinates by the
    /// caller. For the conversion from global cartesian coordinates, the
    /// nested @c shape struct can be used. The point is assumed to be in
    /// the cylinder 2D frame (r * phi, z).
    ///
    /// @tparam is_rad_check whether the radial bound should be checked in this
    /// call.
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

        if constexpr (kRadialCheck) {
            if (std::abs(loc_p[2] - bounds[e_r]) > tol) {
                return false;
            }
        }
        return (bounds[e_n_half_z] - tol <= loc_p[1] and
                loc_p[1] <= bounds[e_p_half_z] + tol);
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
        const scalar_t xy_bound{bounds[e_r] + env};
        return {-xy_bound, -xy_bound, bounds[e_n_half_z] - env,
                xy_bound,  xy_bound,  bounds[e_p_half_z] + env};
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

        if (bounds[e_r] < tol) {
            os << "ERROR: Radius must be in the range (0, numeric_max)"
               << std::endl;
            return false;
        }
        if (bounds[e_n_half_z] >= bounds[e_p_half_z] or
            std::abs(bounds[e_n_half_z] - bounds[e_p_half_z]) < tol) {
            os << "ERROR: Neg. half length must be smaller than pos. half "
                  "length.";
            return false;
        }

        return true;
    }
};

}  // namespace detray
