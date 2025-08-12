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

namespace intersection {

/// Intersection status
enum class status : std::uint_least8_t {
    e_inside = 0u,   ///< Inside the mask within numeric uncertainty
    e_edge = 1u,     ///< Inside the mask's external tolerance band
    e_outside = 2u,  ///< Outside of mask (including all tolerances)
};

}  // namespace intersection

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
    // This is needed for the SIMD implementation, where the boolean type that
    // results from the mask tolerance check has to match the SIMD vector type
    // on which it is used on for a masked assignment of the different status
    // codes
    using status_t = std::conditional_t<algebra::concepts::soa<algebra_t>, T,
                                        intersection::status>;
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

    /// Result of the intersection
    dsimd<algebra_t, status_t> status =
        static_cast<status_t>(intersection::status::e_outside);

    /// Direction of the intersection with respect to the track (true = along,
    /// false = opposite)
    bool_t direction{true};

    /// @returns true if debug information needs to be filled
    static consteval bool is_debug() { return false; }

    /// Set the intersection status according to enum value @param s
    constexpr void set_status(intersection::status s) {
        status = static_cast<status_t>(s);
    }

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

    /// @returns true if any of the intersection results is 'inside'
    DETRAY_HOST_DEVICE
    constexpr bool is_inside() const {
        const dsimd<algebra_t, status_t> comp(
            static_cast<status_t>(intersection::status::e_inside));
        return detail::any_of(this->status == comp);
    }

    /// @returns true if any of the intersection results is 'edge'
    DETRAY_HOST_DEVICE
    constexpr bool is_edge() const {
        const dsimd<algebra_t, status_t> comp(
            static_cast<status_t>(intersection::status::e_edge));
        return detail::any_of(this->status == comp);
    }

    /// @returns true if any of the intersection results is 'inside' or 'edge'
    DETRAY_HOST_DEVICE
    constexpr bool is_probably_inside() const {
        const dsimd<algebra_t, status_t> comp(
            static_cast<status_t>(intersection::status::e_edge));
        return detail::any_of(this->status <= comp);
    }

    /// @returns true if all of the intersection results are 'outside'
    DETRAY_HOST_DEVICE
    constexpr bool is_outside() const {
        const status_t comp(
            static_cast<status_t>(intersection::status::e_outside));
        return detail::all_of(this->status == comp);
    }

    /// Transform to a string for output debugging
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &out_stream,
                                    const intersection2D &is) {
        out_stream << "dist:" << is.path
                   << ", surface: " << is.sf_desc.barcode()
                   << ", type: " << static_cast<int>(is.sf_desc.mask().id())
                   << ", links to vol:" << is.volume_link << ")";

        if (is.is_inside()) {
            out_stream << ", status: inside";
        } else if (is.is_edge()) {
            out_stream << ", status: edge";
        } else {
            out_stream << ", status: outside";
        }
        if constexpr (std::is_scalar_v<bool_t>) {
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

}  // namespace detray
