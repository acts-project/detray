/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// SYCL include(s).
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#include <CL/sycl.hpp>
#elif (IS_SOA)
#include <Vc/Vc>
#endif

// System include(s).
#include <cmath>

namespace detray {

/// Namespace to pick up math functions from
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
namespace math_ns = cl::sycl;
#elif IS_SOA

namespace math_ns {

using std::abs;
using std::asin;
using std::atan;
using std::copysign;
using std::cos;
using std::exp;
using std::fma;
using std::log;
using std::pow;
using std::signbit;
using std::sin;
using std::sqrt;
using std::tan;

/// Vc overloads of common math functions
/// @{
template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) abs(T &&vec) {
    return Vc::abs(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) sqrt(T &&vec) {
    return Vc::sqrt(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) exp(T &&vec) {
    return Vc::exp(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) log(T &&vec) {
    return Vc::log(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) sin(T &&vec) {
    return Vc::sin(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) asin(T &&vec) {
    return Vc::asin(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) cos(T &&vec) {
    return Vc::cos(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) tan(T &&vec) {
    // It seems there is no dedicated @c Vc::tan function ?
    return Vc::sin(std::forward<T>(vec)) / Vc::cos(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) atan(T &&vec) {
    return Vc::atan(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) copysign(T &&mag, T &&sgn) {
    return Vc::copysign(std::forward<T>(mag), std::forward<T>(sgn));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) signbit(T &&vec) {
    return Vc::isnegative(std::forward<T>(vec));
}

template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_vector<T>::value, bool> = true>
inline decltype(auto) fma(T &&x, T &&y, T &&z) {
    return Vc::fma(std::forward<T>(x), std::forward<T>(y), std::forward<T>(z));
}
/// @}

}  // namespace math_ns

#else
namespace math_ns = std;
#endif  // SYCL

namespace detail {

/// boolean utilities
/// @{
inline constexpr bool any_of(bool b) {
    // Do nothing
    return b;
}
#if IS_SOA
template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_mask<T>::value, bool> = true>
inline bool any_of(T &&mask) {
    return Vc::any_of(std::forward<T>(mask));
}
#endif
/// @}

using math_ns::copysign;
using math_ns::signbit;

/// Composes a floating point value with the magnitude of @param mag and the
/// sign of @param sgn
/*template <typename scalar_t>
DETRAY_HOST_DEVICE inline scalar_t copysign(scalar_t mag, scalar_t sgn) {
#if defined(__CUDACC__)
    if constexpr (std::is_same_v<scalar_t, float>) {
        return copysignf(mag, sgn);
    } else {
        return copysign(mag, sgn);
    }
#elif !defined(__CUDACC__)
    return math_ns::copysign(mag, sgn);
#endif
}

/// Gets the signbit from a variable
template <typename scalar_t>
DETRAY_HOST_DEVICE inline bool signbit(scalar_t arg) {
#if defined(__CUDACC__)
    return signbit(arg);
#elif !defined(__CUDACC__)
    return math_ns::signbit(arg);
#endif
}*/

}  // namespace detail

}  // namespace detray
