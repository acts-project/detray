/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
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
#elif DETRAY_ALGEBRA_VC_AOS
#include "detray/plugins/algebra/vc_aos_definitions.hpp"
#elif DETRAY_ALGEBRA_VC_SOA
#include "detray/plugins/algebra/vc_soa_definitions.hpp"
#else
#error "No algebra plugin selected! Please link to one of the algebra plugins."
#endif

// Algebra-plugins include(s)
#include "algebra/utils/print.hpp"

namespace detray {

// Pull in the print operator definitions for the algebra types
using algebra::operator<<;

template <typename A>
using dvalue = typename algebra::get_value_t<A>;

template <typename A>
using dbool = typename algebra::get_boolean_t<A>;

template <typename A, typename T>
using dsimd = algebra::get_simd_t<A, T>;

template <typename A>
using dsize_type = algebra::get_size_t<A>;

template <typename A>
using dscalar = algebra::get_scalar_t<A>;

template <typename A>
using dpoint2D = algebra::get_point2D_t<A>;

template <typename A>
using dpoint3D = algebra::get_point3D_t<A>;

template <typename A>
using dvector3D = algebra::get_vector3D_t<A>;

template <typename A>
using dtransform3D = algebra::get_transform3D_t<A>;

template <typename A, std::size_t R, std::size_t C>
using dmatrix = algebra::get_matrix_t<A, R, C>;

namespace detail {

using namespace ::algebra::boolean;

}  // namespace detail

}  // namespace detray
