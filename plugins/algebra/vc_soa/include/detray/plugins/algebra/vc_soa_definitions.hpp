/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/vc_soa.hpp"

#define __plugin algebra::vc_soa
#define ALGEBRA_PLUGIN vc_soa
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
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = Vc::Vector<T>;

    using boolean = Vc::Mask<V>;
    using scalar = simd<value_type>;
    using transform3D = algebra::vc_soa::transform3<value_type>;
    using point3D = algebra::vc_soa::point3<value_type>;
    using vector3D = algebra::vc_soa::vector3<value_type>;
};
/// @}

// Define namespace(s)
namespace getter = algebra::getter;
namespace vector = algebra::vector;
namespace matrix = algebra::matrix;

}  // namespace detray
