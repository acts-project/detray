/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if DETRAY_ALGEBRA_ARRAY
#include "detray/plugins/algebra/array_definitions.hpp"
#elif DETRAY_ALGEBRA_EIGEN
#include "detray/plugins/algebra/eigen_definitions.hpp"
#elif DETRAY_ALGEBRA_SMATRIX
#include "detray/plugins/algebra/smatrix_definitions.hpp"
#elif DETRAY_ALGEBRA_VC
#include "detray/plugins/algebra/vc_array_definitions.hpp"
#elif DETRAY_ALGEBRA_VC_SOA
#include "detray/plugins/algebra/vc_soa_definitions.hpp"
#else
#error "No algebra plugin selected! Please link to one of the algebra plugins."
#endif
