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

// System include(s)
#include <type_traits>

namespace detray {

namespace detail {
/// The detray scalar/boolean types (can be SIMD)
/// @{
template <typename T, typename = void>
struct get_scalar {};

template <typename T>
struct get_scalar<T, std::enable_if_t<std::is_arithmetic_v<T>, void>> {
    using scalar = T;
    using boolean = bool;
};

template <typename T>
struct get_scalar<
    T, std::enable_if_t<!std::is_same_v<typename T::scalar, void>, void>> {
    using scalar = typename T::scalar;
    using boolean = typename T::boolean;
};
/// @}

/// The detray algebra types (can be SIMD)
/// @{
template <typename T, typename = void>
struct get_algebra {};

template <typename T>
struct get_algebra<T, std::enable_if_t<std::is_arithmetic_v<T>, void>> {
    using point3D = std::array<T, 3>;
    using vector3D = std::array<T, 3>;
    using transform3D = std::array<T, 9>;
};

template <typename T>
struct get_algebra<
    T, std::enable_if_t<!std::is_same_v<typename T::point3D, void>, void>> {
    using point2D = typename T::point2D;
    using point3D = typename T::point3D;
    using vector3D = typename T::vector3D;
    using transform3D = typename T::transform3D;
};
/// @}
}  // namespace detail

template <template <typename> class A, typename T>
using dsimd = typename A<float>::simd<T>;

template <typename A = detray::scalar>
using dscalar = typename detail::get_scalar<A>::scalar;

template <typename A>
using dpoint2D = typename detail::get_algebra<A>::point2D;

template <typename A>
using dpoint3D = typename detail::get_algebra<A>::point3D;

template <typename A>
using dvector3D = typename detail::get_algebra<A>::vector3D;

template <typename A>
using dtransform3D = typename detail::get_algebra<A>::transform3D;

template <typename A>
using dbool = typename detail::get_scalar<A>::boolean;

}  // namespace detray
