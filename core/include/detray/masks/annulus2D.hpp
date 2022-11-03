/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/polar2.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_bounds.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <string>
#include <tuple>

namespace detray {

/// @brief Geometrical shape of a stereo annulus that is used for the itk
/// strip endcaps.
///
/// @tparam intersector_t defines how to intersect the underlying surface
///         geometry
///
/// The stereo annulus is defined in two different(!) polar coordinate systems
/// that differ by an origin shift. The boundaries are the inner and outer
/// radius (bounds[0] and bounds[1]) in the polar coordinate system of an
/// endcap strip disc (called disc system in the following) that is centered on
/// the beam axis, as well as the two phi boundaries that are defined in the
/// system that is shifted into the focal point of the strips on a given module
/// (called focal system in the following). The mask phi boundaries (bounds[2]
/// and bounds[3]) are defined relative to the average phi position of the
/// strips in the focal system. Using a conversion between the two coordinate
/// systems, these boundaries can be checked with a tolerance in r (tol[0]), as
/// well as in phi (tol[1]).
/// Due to the focal polar coordinate system of the strips needing a
/// different origin from the discs polar system, three additional conversion
/// parameters are included (bounds[4], bounds[5], bounds[6]).
/// The first two are the origin shift in x and y respectively, while
/// bounds[6] is the average Phi angle mentioned above.
template <template <typename> class intersector_t = plane_intersector>
class annulus2D {
    public:
    /// The name for this shape
    inline static const std::string name = "(stereo) annulus2D";

    /// Names for the mask boundary values
    enum boundaries : std::size_t {
        e_min_r = 0,
        e_max_r = 1,
        e_min_phi_rel = 2,
        e_max_phi_rel = 3,
        e_shift_x = 4,
        e_shift_y = 5,
        e_average_phi = 6,
        e_size = 7,
    };

    /// Local coordinate frame (both for disc and focal system ?)
    template <typename algebra_t>
    using local_frame_type = polar2<algebra_t>;
    /// Local point type (2D)
    template <typename algebra_t>
    using loc_point_type = typename local_frame_type<algebra_t>::point2;

    /// Measurement frame. @todo switch to focal system?
    template <typename algebra_t>
    using measurement_frame_type = polar2<algebra_t>;
    /// Local measurement point (2D)
    template <typename algebra_t>
    using measurement_point_type =
        typename measurement_frame_type<algebra_t>::point2;

    /// Underlying surface geometry: planar
    template <typename algebra_t>
    using intersector_type = intersector_t<algebra_t>;

    /// Behaviour of the two local axes (linear in r, circular in phi)
    template <
        n_axis::bounds e_s = n_axis::bounds::e_closed,
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_r;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_phi;

        using types = std::tuple<n_axis::bounds_t<e_s, axis_loc0>,
                                 n_axis::circular<axis_loc1>>;

        /// Local coordinate frame (both for disc and focal system ?)
        template <typename algebra_t>
        using coordinate_type = local_frame_type<algebra_t>;

        template <typename C, typename S>
        using binning = std::tuple<binning_loc0<C, S>, binning_loc1<C, S>>;
    };

    /// Given a local polar point given in the disc frame, @returns the
    /// correponding point in the focal frame
    /*DETRAY_HOST_DEVICE
    template<typename scalar_t>
    static inline constexpr loc_point_type<scalar_t> to_focal_frame(
        const loc_point_type<scalar_t> & pc_mod_point) {
        return {};
    }

    /// Given a local polar point given in the focal frame, @returns the
    /// correponding point in the disc frame
    DETRAY_HOST_DEVICE
    template<typename scalar_t>
    static inline constexpr loc_point_type<scalar_t> to_disc_frame(
        const loc_point_type<scalar_t> & pc_strp_point) {
        return {};
    }*/

    /// @brief Check boundary values for a local point.
    ///
    /// @note the point is expected to be given in local coordinates by the
    /// caller. For the conversion from global cartesian coordinates, the
    /// nested @c shape struct can be used. For this mask, the local point must
    /// be given in the focal (strips) system.
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

        // The two quantities to check: r^2 in disc system, phi in focal system:

        // Rotate by avr phi into the focal system
        const scalar_t phi_strp{loc_p[1] - bounds[e_average_phi]};

        // Check phi boundaries, which are well def. in focal frame
        if ((phi_strp < bounds[e_min_phi_rel] - tol) or
            (phi_strp > bounds[e_max_phi_rel] + tol)) {
            return false;
        }

        // Now go to module frame to check r boundaries. Use the origin
        // shift in polar coordinates for that
        // TODO: Put shift in r-phi into the bounds?
        const point_t shift_xy = {scalar_t{-1} * bounds[e_shift_x],
                                  scalar_t{-1} * bounds[e_shift_y]};
        const scalar_t shift_r{getter::perp(shift_xy)};
        const scalar_t shift_phi{getter::phi(shift_xy)};

        const scalar_t r_mod2{shift_r * shift_r + loc_p[0] * loc_p[0] +
                              scalar_t{2} * shift_r * loc_p[0] *
                                  std::cos(phi_strp - shift_phi)};

        // Apply tolerances
        const scalar_t minR_tol{bounds[e_min_r] - tol};
        const scalar_t maxR_tol{bounds[e_max_r] + tol};

        return ((r_mod2 >= minR_tol * minR_tol) and
                (r_mod2 <= maxR_tol * maxR_tol));
    }
};

}  // namespace detray
