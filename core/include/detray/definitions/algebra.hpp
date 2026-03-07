/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if DETRAY_ALGEBRA_ARRAY
#include "algebra/array.hpp"
#elif DETRAY_ALGEBRA_EIGEN
#include "algebra/eigen.hpp"
#elif DETRAY_ALGEBRA_FASTOR
#include "algebra/fastor.hpp"
#elif DETRAY_ALGEBRA_SMATRIX
#include "algebra/smatrix.hpp"
#elif DETRAY_ALGEBRA_VC_AOS
#include "algebra/vc_aos.hpp"
#elif DETRAY_ALGEBRA_VC_SOA
#include "algebra/vc_soa.hpp"
#else
#error "No algebra plugin selected! Please link to one of the algebra plugins."
#endif

// Algebra-plugins include(s)
#include "detray/algebra/common/boolean.hpp"
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/algebra/utils/casts.hpp"
#include "detray/algebra/utils/print.hpp"

namespace detray {

template <typename A>
using dvalue = get_value_t<A>;

template <typename A>
using dbool = get_boolean_t<A>;

template <typename A, typename T>
using dsimd = get_simd_t<A, T>;

template <typename A>
using dsize_type = get_size_t<A>;

template <typename A>
using dscalar = get_scalar_t<A>;

template <typename A>
using dpoint2D = get_point2D_t<A>;

template <typename A>
using dpoint3D = get_point3D_t<A>;

template <typename A>
using dvector3D = get_vector3D_t<A>;

template <typename A>
using dtransform3D = get_transform3D_t<A>;

template <typename A, std::size_t R, std::size_t C>
using dmatrix = get_matrix_t<A, R, C>;

}  // namespace detray
