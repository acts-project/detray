/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray core include(s).
#include "detray/definitions/algebra.hpp"

// Detray detector include(s)
#include "detray/detectors/default_metadata.hpp"
#include "detray/detectors/toy_metadata.hpp"

namespace detray::benchmarks {

// Select algebra-plugin to compile the test with
#if DETRAY_ALGEBRA_ARRAY
using algebra = detray::array<DETRAY_CUSTOM_SCALARTYPE>;

#elif DETRAY_ALGEBRA_EIGEN
using algebra = detray::eigen<DETRAY_CUSTOM_SCALARTYPE>;

#elif DETRAY_ALGEBRA_SMATRIX
using algebra = detray::smatrix<DETRAY_CUSTOM_SCALARTYPE>;

#elif DETRAY_ALGEBRA_VC_AOS
using algebra = detray::vc_aos<DETRAY_CUSTOM_SCALARTYPE>;
#endif

// Test algebra types
using scalar = dscalar<algebra>;
using point2 = dpoint2D<algebra>;
using point3 = dpoint3D<algebra>;
using vector3 = dvector3D<algebra>;
using transform3 = dtransform3D<algebra>;
template <std::size_t ROWS, std::size_t COLS>
using matrix = dmatrix<algebra, ROWS, COLS>;

// Test detector types
using default_metadata = detray::default_metadata<algebra>;
using toy_metadata = detray::toy_metadata<algebra>;

}  // namespace detray::benchmarks
