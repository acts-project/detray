/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/polar2.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_bounds.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <string>

namespace detray {

/// @brief Geometrical shape of a stereo annulus that is used for the itk
/// strip endcaps.
///
/// @tparam intersector_t defines how to intersect the underlying surface
///         geometry
/// @tparam kMeasDim defines the dimension of the measurement
/// @tparam kNormalOrder true if the index for measurement parameter follows
/// the local coordinate system
///
/// The stereo annulus is defined in two different(!) polar coordinate systems
/// that differ by an origin shift. The boundaries are the inner and outer
/// radius (bounds[0] and bounds[1]) in the polar coordinate system of an
/// endcap disc (called beam system in the following) that is centered on
/// the beam axis, as well as the two phi boundaries that are defined in the
/// system that is shifted into the focal point of the strips on a given sensor
/// (called focal system in the following). Note, that the local coordinate
/// system of the annulus surface is the same as the shifted disc (focal)
/// system! The mask phi boundaries (bounds[2] and bounds[3]) are defined
/// relative to the average phi position ( bounds[6]) of the strips in the
/// focal system.
/// Due to the focal polar coordinate system of the strips needing a different
/// origin from the beam polar system, two additional conversion parameters are
/// included (bounds[4], bounds[5]). These are the origin shift in x and y
/// respectively.
template <template <typename> class intersector_t = plane_intersector,
          unsigned int kMeasDim = 1u, bool kNormalOrder = false>
class annulus2D {
    public:
    /// The name for this shape
    inline static const std::string name = "(stereo) annulus2D";

    /// The measurement dimension
    inline static constexpr const unsigned int meas_dim{kMeasDim};

    /// Normal ordering
    inline static constexpr const bool normal_order{kNormalOrder};

    // Measurement dimension check
    static_assert(meas_dim == 1u || meas_dim == 2u,
                  "Only 1D or 2D measurement is allowed");

    /// Names for the mask boundary values
    enum boundaries : unsigned int {
        e_min_r = 0u,
        e_max_r = 1u,
        e_min_phi_rel = 2u,
        e_max_phi_rel = 3u,
        e_average_phi = 4u,
        e_shift_x = 5u,
        e_shift_y = 6u,
        e_size = 7u,
    };

    /// Local coordinate frame ( focal system )
    template <typename algebra_t>
    using local_frame_type = polar2<algebra_t>;

    /// Underlying surface geometry: planar
    template <typename intersection_t>
    using intersector_type = intersector_t<intersection_t>;

    /// Behaviour of the two local axes (linear in r, circular in phi)
    template <
        n_axis::bounds e_s = n_axis::bounds::e_closed,
        template <typename, typename> class binning_loc0 = n_axis::regular,
        template <typename, typename> class binning_loc1 = n_axis::regular>
    struct axes {
        static constexpr n_axis::label axis_loc0 = n_axis::label::e_r;
        static constexpr n_axis::label axis_loc1 = n_axis::label::e_phi;
        static constexpr std::size_t dim{2u};

        using types = dtuple<n_axis::bounds_t<e_s, axis_loc0>,
                             n_axis::circular<axis_loc1>>;

        /// Local coordinate frame (both for disc and focal system ?)
        template <typename algebra_t>
        using coordinate_type = local_frame_type<algebra_t>;

        template <typename C, typename S>
        using binning = dtuple<binning_loc0<C, S>, binning_loc1<C, S>>;
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

    /// @returns the stereo angle calculated from the mask @param bounds .
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE std::array<scalar_t, 8> stereo_angle(
        const bounds_t<scalar_t, kDIM>& bounds) const {
        // Half stereo angle (phi_s / 2) (y points in the long strip direction)
        return 2.f * math_ns::atan(bounds[e_shift_y] / bounds[e_shift_x]);
    }

    /// @brief Check boundary values for a local point.
    ///
    /// @note the point is expected to be given in local coordinates by the
    /// caller. For the annulus shape, the local coordinate system of the
    /// strips is used.
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
        const bounds_t<scalar_t, kDIM>& bounds, const point_t& loc_p,
        const scalar_t tol = std::numeric_limits<scalar_t>::epsilon()) const {

        // The two quantities to check: r^2 in beam system, phi in focal system:

        // Rotate by avr phi in the focal system (this is usually zero)
        const scalar_t phi_strp{loc_p[1] - bounds[e_average_phi]};

        // Check phi boundaries, which are well def. in focal frame
        if ((phi_strp < bounds[e_min_phi_rel] - tol) or
            (phi_strp > bounds[e_max_phi_rel] + tol)) {
            return false;
        }

        // Now go to beam frame to check r boundaries. Use the origin
        // shift in polar coordinates for that
        // TODO: Put shift in r-phi into the bounds?
        const point_t shift_xy = {-bounds[e_shift_x], -bounds[e_shift_y], 0.f};
        const scalar_t shift_r{getter::perp(shift_xy)};
        const scalar_t shift_phi{getter::phi(shift_xy)};

        const scalar_t r_mod2{shift_r * shift_r + loc_p[0] * loc_p[0] +
                              2.f * shift_r * loc_p[0] *
                                  math_ns::cos(phi_strp - shift_phi)};

        // Apply tolerances as squares: 0 <= a, 0 <= b: a^2 <= b^2 <=> a <= b
        const scalar_t minR_tol{bounds[e_min_r] - tol};
        const scalar_t maxR_tol{bounds[e_max_r] + tol};

        assert(minR_tol >= 0.f);

        return ((r_mod2 >= minR_tol * minR_tol) and
                (r_mod2 <= maxR_tol * maxR_tol));
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
    // @TODO: this is a terrible approximation: restrict to annulus corners
    template <typename algebra_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE std::array<scalar_t, 6> local_min_bounds(
        const bounds_t<scalar_t, kDIM>& bounds,
        const scalar_t env = std::numeric_limits<scalar_t>::epsilon()) const {

        using point_t = typename algebra_t::point2;

        assert(env > 0.f);

        const auto c_pos = corners(bounds);

        const scalar_t o_x{bounds[e_shift_x]};
        const scalar_t o_y{bounds[e_shift_y]};

        // Corner points 'b' and 'c' in local cartesian beam system
        const point_t b{c_pos[4] * math_ns::cos(c_pos[5]) - o_x,
                        c_pos[4] * math_ns::sin(c_pos[5]) - o_y};
        const point_t c{c_pos[6] * math_ns::cos(c_pos[7]) - o_x,
                        c_pos[6] * math_ns::sin(c_pos[7]) - o_y};

        // bisector = 0.5 * (c + b). Scale to the length of the full circle to
        // get the outermost point
        const point_t t = bounds[e_max_r] * vector::normalize(c + b);

        // Find the min/max positions in x and y
        std::array<scalar_t, 5> x_pos{
            c_pos[2] * math_ns::cos(c_pos[3]) - o_x, b[0], c[0],
            c_pos[0] * math_ns::cos(c_pos[1]) - o_x, t[0]};
        std::array<scalar_t, 5> y_pos{
            c_pos[2] * math_ns::sin(c_pos[3]) - o_y, b[1], c[1],
            c_pos[0] * math_ns::sin(c_pos[1]) - o_y, t[1]};

        constexpr scalar_t inf{std::numeric_limits<scalar_t>::infinity()};
        scalar_t min_x{inf}, min_y{inf}, max_x{-inf}, max_y{-inf};
        for (unsigned int i{0u}; i < 5u; ++i) {
            min_x = x_pos[i] < min_x ? x_pos[i] : min_x;
            max_x = x_pos[i] > max_x ? x_pos[i] : max_x;
            min_y = y_pos[i] < min_y ? y_pos[i] : min_y;
            max_y = y_pos[i] > max_y ? y_pos[i] : max_y;
        }

        return {min_x - env, min_y - env, -env, max_x + env, max_y + env, env};
    }

    /// @brief Stereo annulus corners in polar strip system.
    ///
    /// @param bounds the boundary values for the stereo annulus
    ///
    /// @note see caluclation of strip lengths in
    ///       https://cds.cern.ch/record/1514636?ln=en p10-13
    ///
    /// @returns an array of coordinates that contains the lower point (first
    /// four values) and the upper point (latter four values).
    template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST_DEVICE std::array<scalar_t, 8> corners(
        const bounds_t<scalar_t, kDIM>& bounds) const {

        // Calculate the r-coordinate of a point in the strip system from the
        // circle arc radius (e.g. min_r) and the phi position in the strip
        // system (e.g. for the corners these are average_phi + min_phi_rel and
        // average_phi + max_phi_rel).
        auto get_strips_pc_r = [&bounds](const scalar_t R,
                                         const scalar_t phi) -> scalar_t {
            // f: Shift distance between beamline and focal system origin
            const scalar_t f2{bounds[e_shift_x] * bounds[e_shift_x] +
                              bounds[e_shift_y] * bounds[e_shift_y]};

            // f * sin(phi_s / 2 + phi) using: f_y / f = sin(phi_s / 2) and
            // sin(a + b) = sin(a)*cos(b) + cos(a)*sin(b)
            const scalar_t f_sin_phi{bounds[e_shift_x] * math_ns::cos(phi) +
                                     bounds[e_shift_y] * math_ns::sin(phi)};

            return f_sin_phi +
                   math_ns::sqrt(f_sin_phi * f_sin_phi - f2 + R * R);
        };

        // Calculate the polar coordinates for the corners
        const scalar_t min_phi{bounds[e_average_phi] + bounds[e_min_phi_rel]};
        const scalar_t max_phi{bounds[e_average_phi] + bounds[e_max_phi_rel]};
        std::array<scalar_t, 8> corner_pos;
        // bottom left: min_r, min_phi_rel
        corner_pos[0] = get_strips_pc_r(bounds[e_min_r], min_phi);
        corner_pos[1] = min_phi;
        // bottom right: min_r, max_phi_rel
        corner_pos[2] = get_strips_pc_r(bounds[e_min_r], max_phi);
        corner_pos[3] = max_phi;
        // top right: max_r, max_phi_rel
        corner_pos[4] = get_strips_pc_r(bounds[e_max_r], max_phi);
        corner_pos[5] = max_phi;
        // top left: max_r, min_phi_rel
        corner_pos[6] = get_strips_pc_r(bounds[e_max_r], min_phi);
        corner_pos[7] = min_phi;

        return corner_pos;
    }

    /// @brief Calculates the coordinates of the vertices.
    ///
    /// @param bounds the boundary values for this shape.
    ///
    /// @returns a container of vertices. If the shape contains
    /// no vertices an empty container will be returned.
    template <typename point3_container_t,
              template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST inline point3_container_t local_vertices(
        const bounds_t<scalar_t, kDIM>&) const {
        return {};
    }

    /// @brief Finds the closest point lying on the surface to the given point.
    ///
    /// @param bounds the boundary values for this shape.
    /// @param loc_p the point in the local coordinate system.
    ///
    /// @returns the closest point lying on the surface in the local_coordinate system.
        template <template <typename, std::size_t> class bounds_t,
              typename scalar_t, std::size_t kDIM, typename point_t,
              typename std::enable_if_t<kDIM == e_size, bool> = true>
    DETRAY_HOST inline point_t closest_surface_point(
        const bounds_t<scalar_t, kDIM>& bounds, const point_t& loc_p) const {
            return point_t{};
    }
};

}  // namespace detray
