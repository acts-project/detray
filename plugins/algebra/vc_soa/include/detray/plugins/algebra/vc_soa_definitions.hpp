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

/// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
using transform3D = algebra::vc_soa::transform3<scalar>;
using point3D = algebra::vc_soa::point3<scalar>;
using vector3D = algebra::vc_soa::vector3<scalar>;
/// @}

/// 

// Define namespace(s)
namespace getter = algebra::getter;
namespace vector = algebra::vector;
namespace matrix = algebra::matrix;

}  // namespace detray
