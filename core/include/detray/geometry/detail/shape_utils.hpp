/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"

namespace detray::detail {

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

}  // namespace detray::detail
