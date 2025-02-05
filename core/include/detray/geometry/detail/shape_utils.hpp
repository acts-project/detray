// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"

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
