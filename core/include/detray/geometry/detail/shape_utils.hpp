/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {

namespace detail {

/// Generate phi tolerance from distance tolerance
///
/// @param tol is the distance tolerance in mm
/// @param radius is the radius of the shape
///
/// @return the opening angle of a chord the size of tol (= 2*arcsin(c/(2r)))
/// using a small angle approximation
template <concepts::scalar scalar_t>
constexpr scalar_t phi_tolerance(scalar_t tol, scalar_t radius) {
    return radius > 0.f ? tol / radius : tol;
}
// Result of an 'inside' check
template <typename bool_t>
using boundary_check_result = darray<bool_t, 2>;

}  // namespace detail

/// Address different types of check results
enum class check_type : std::uint_least8_t {
    e_inside = 0u,
    e_with_edge = 1u,
};

/// Get the result of a check with edge tolerance
template <check_type C, typename bool_t>
DETRAY_HOST_DEVICE constexpr bool_t get(
    const detail::boundary_check_result<bool_t> result) {
    if constexpr (C == check_type::e_inside) {
        return result[0];
    } else if constexpr (C == check_type::e_with_edge) {
        return result[1];
    } else {
        // Broadcast
        return bool_t(false);
    }
}

/// Get the result of a simple tolerance check
template <check_type C, typename bool_t>
DETRAY_HOST_DEVICE constexpr bool_t get(const bool_t& result) {
    return result;
}

}  // namespace detray
