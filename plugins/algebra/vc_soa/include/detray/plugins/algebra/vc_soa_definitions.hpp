/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/math/vc_soa.hpp"
#include "algebra/storage/vc_soa.hpp"

#define IS_SOA 1

namespace detray {

using algebra::storage::operator*;
using algebra::storage::operator/;
using algebra::storage::operator-;
using algebra::storage::operator+;

using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct vc_soa {
    /// Define scalar precision
    using value_type = V;

    template <typename T>
    using simd = Vc::Vector<T>;

    using boolean = Vc::Mask<V>;

    /// Linear Algebra type definitions
    /// @{
    using scalar = simd<value_type>;
    using transform3D =
        algebra::vc_soa::math::transform3<algebra::vc_soa::storage_type,
                                          value_type>;
    using point2D = algebra::vc_soa::point2<value_type>;
    using point3D = algebra::vc_soa::point3<value_type>;
    using vector3D = algebra::vc_soa::vector3<value_type>;
    /// @}
};
/// @}

namespace vector {

using algebra::vc_soa::math::cross;
using algebra::vc_soa::math::dot;
using algebra::vc_soa::math::normalize;

}  // namespace vector

namespace getter {

using algebra::vc_soa::math::eta;
using algebra::vc_soa::math::norm;
using algebra::vc_soa::math::perp;
using algebra::vc_soa::math::phi;
using algebra::vc_soa::math::theta;

}  // namespace getter

}  // namespace detray
