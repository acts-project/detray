/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <cstdint>
#include <limits>
#include <ostream>

namespace detray {

template <typename surface_descr_t, concepts::algebra algebra_t,
          bool debug = false>
struct intersection2D {};

/// @brief This class holds the intersection information.
///
/// @tparam surface_descr_t is the type of surface descriptor
template <typename surface_descr_t, concepts::algebra algebra_t>
struct intersection2D<surface_descr_t, algebra_t, false> {

    using algebra_type = algebra_t;
    using T = dvalue<algebra_t>;
    using bool_t = dbool<algebra_t>;
    using scalar_type = dscalar<algebra_t>;
    using point2_type = dpoint2D<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    using nav_link_t = typename surface_descr_t::navigation_link;

    /// Descriptor of the surface this intersection belongs to
    surface_descr_t sf_desc{};

    /// Distance between track and candidate
    scalar_type path = detail::invalid_value<T>();

    /// Navigation information (next volume to go to)
    nav_link_t volume_link{detail::invalid_value<nav_link_t>()};

    /// Result of the intersection (true = inside, false = outside)
    bool_t status{false};

    /// Direction of the intersection with respect to the track (true = along,
    /// false = opposite)
    bool_t direction{true};

    /// @returns true if debug information needs to be filled
    static consteval bool is_debug() { return false; }

    /// @note: Three way comparison cannot be used easily with SoA boolean masks
    /// @{
    /// @param rhs is the right hand side intersection for comparison
    DETRAY_HOST_DEVICE
    friend constexpr bool_t operator<(const intersection2D &lhs,
                                      const intersection2D &rhs) noexcept {
        return (math::fabs(lhs.path) < math::fabs(rhs.path));
    }

    /// @param rhs is the right hand side intersection for comparison
    DETRAY_HOST_DEVICE
    friend constexpr bool_t operator<=(const intersection2D &lhs,
                                       const intersection2D &rhs) noexcept {
        return (math::fabs(lhs.path) <= math::fabs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    friend constexpr bool_t operator>(const intersection2D &lhs,
                                      const intersection2D &rhs) noexcept {
        return (math::fabs(lhs.path) > math::fabs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    friend constexpr bool_t operator>=(const intersection2D &lhs,
                                       const intersection2D &rhs) noexcept {
        return (math::fabs(lhs.path) > math::fabs(rhs.path));
    }
    /// @}

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    friend constexpr bool_t operator==(const intersection2D &lhs,
                                       const intersection2D &rhs) noexcept {
        return math::fabs(lhs.path - rhs.path) <
               std::numeric_limits<float>::epsilon();
    }

    DETRAY_HOST_DEVICE
    constexpr bool_t is_inside() const { return detail::any_of(this->status); }

    /// Transform to a string for output debugging
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &out_stream,
                                    const intersection2D &is) {
        out_stream << "dist:" << is.path
                   << ", surface: " << is.sf_desc.barcode()
                   << ", type: " << static_cast<int>(is.sf_desc.mask().id())
                   << ", links to vol:" << is.volume_link << ")";

        if constexpr (std::is_scalar_v<bool_t>) {
            out_stream << (is.status ? ", status: inside"
                                     : ", status: outside");
            out_stream << (is.direction ? ", direction: along"
                                        : ", direction: opposite");
        } else {
            out_stream << ", status: " << is.status;
            out_stream << ", direction: " << is.direction;
        }

        return out_stream;
    }
};

/// @brief This class holds the intersection information with additional debug
///        information.
///
/// @tparam surface_descr_t is the type of surface descriptor
template <typename surface_descr_t, concepts::algebra algebra_t>
struct intersection2D<surface_descr_t, algebra_t, true>
    : public intersection2D<surface_descr_t, algebra_t, false> {

    using T = dvalue<algebra_t>;
    using algebra_type = algebra_t;
    using point3_type = dpoint3D<algebra_t>;

    /// @returns true if debug information needs to be filled
    static consteval bool is_debug() { return true; }

    /// Local position of the intersection on the surface
    point3_type local{detail::invalid_value<T>(), detail::invalid_value<T>(),
                      detail::invalid_value<T>()};

    /// Transform to a string for output debugging
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &out_stream,
                                    const intersection2D &is) {
        using base_t = intersection2D<surface_descr_t, algebra_t, false>;

        out_stream << static_cast<base_t>(is);
        out_stream << ", loc [" << is.local[0] << ", " << is.local[1] << ", "
                   << is.local[2] << "], ";

        return out_stream;
    }
};

/// @brief Fill an intersection with the result of the intersection alg.
///
/// @param [in] traj the test trajectory that intersects the surface
/// @param [out] sfi the surface intersection
/// @param [in] s path length to the intersection point
/// @param [in] mask the mask of the surface
/// @param [in] trf the transform of the surface
/// @param [in] mask_tol_scalor scale factor for adaptive mask tolerance
/// @param [in] mask_tolerance minimal and maximal mask tolerance
template <typename trajectory_t, typename intersection_t,
          concepts::point point_t, concepts::scalar scalar_t,
          typename surface_descr_t, typename mask_t,
          concepts::transform3D transform3_t>
DETRAY_HOST_DEVICE constexpr void build_intersection(
    const trajectory_t &traj, intersection_t &is, const point_t &pos,
    const scalar_t s, const surface_descr_t sf_desc, const mask_t &mask,
    const transform3_t &trf, const darray<scalar_t, 2> &mask_tolerance,
    const scalar_t mask_tol_scalor = 0.f, const scalar_t overstep_tol = 0.f) {

    using algebra_t = typename intersection_t::algebra_type;

    // Build intersection struct from test trajectory, if the distance is valid
    if (detray::detail::none_of(s >= overstep_tol)) {
        // Not a valid intersection
        is.status = decltype(is.status)(false);
        return;
    }

    // Save local position for debug evaluation
    if constexpr (intersection_t::is_debug()) {
        // Global position on the surface
        dpoint3D<algebra_t> glob_pos;
        if constexpr (concepts::soa<algebra_t>) {
            // The trajectory is given in AoS layout. Translate...
            const auto &origin = traj.pos();
            const auto &dir = traj.dir();

            // Broadcast
            const dvector3D<algebra_t> ro{origin[0], origin[1], origin[2]};
            const dvector3D<algebra_t> rd{dir[0], dir[1], dir[2]};

            glob_pos = ro + s * rd;
        } else {
            // Works for any parameterized trajectory
            glob_pos = traj.pos(s);
        }
        is.local = mask_t::to_local_frame3D(trf, glob_pos, traj.dir(s));
    }

    // Tol.: scale with distance of surface to account for track bending
    const scalar_t base_tol = math::max(
        mask_tolerance[0],
        math::min(mask_tolerance[1], mask_tol_scalor * math::fabs(s)));

    // Intersector provides specialized local point
    if constexpr (std::same_as<point_t, dpoint2D<algebra_t>>) {
        is.status = mask.is_inside(pos, base_tol);
    } else {
        // Otherwise, let the shape transform the point to local
        is.status = mask.is_inside(trf, pos, base_tol);
    }
    // Mask out solutions that don't meet the overstepping tolerance (SoA)
    is.status &= (is.path >= overstep_tol);

    is.path = s;
    is.sf_desc = sf_desc;
    is.direction = !detail::signbit(is.path);
    is.volume_link = mask.volume_link();
}

}  // namespace detray
