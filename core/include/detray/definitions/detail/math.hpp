/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <concepts>

// SYCL include(s).
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#include <CL/sycl.hpp>
#endif

// System include(s).
#include <cmath>

namespace detray {

/// Namespace to pick up math functions from
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
namespace math = cl::sycl;
#elif IS_SOA

namespace math {

using std::abs;
using std::asin;
using std::atan;
using std::copysign;
using std::cos;
using std::exp;
using std::fabs;
using std::fma;
using std::log;
using std::max;
using std::min;
using std::pow;
using std::signbit;
using std::sin;
using std::sqrt;
using std::tan;

/// Vc overloads of common math functions
/// @{
template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) abs(
    T &&vec) {
    return Vc::abs(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) fabs(
    T &&vec) {
    return Vc::abs(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) sqrt(
    T &&vec) {
    return Vc::sqrt(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) exp(
    T &&vec) {
    return Vc::exp(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) log(
    T &&vec) {
    return Vc::log(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) sin(
    T &&vec) {
    return Vc::sin(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) asin(
    T &&vec) {
    return Vc::asin(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) cos(
    T &&vec) {
    return Vc::cos(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) tan(
    T &&vec) {
    // It seems there is no dedicated @c Vc::tan function ?
    return Vc::sin(std::forward<T>(vec)) / Vc::cos(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) atan(
    T &&vec) {
    return Vc::atan(std::forward<T>(vec));
}

template <typename T, typename S>
requires Vc::Traits::is_simd_vector<T>::value
    &&Vc::Traits::is_simd_vector<S>::value inline decltype(auto)
    copysign(T &&mag, S &&sgn) {
    return Vc::copysign(std::forward<T>(mag), std::forward<S>(sgn));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) min(
    T &&vec) {
    return Vc::min(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) max(
    T &&vec) {
    return Vc::max(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) signbit(
    T &&vec) {
    return Vc::isnegative(std::forward<T>(vec));
}

template <typename T>
requires Vc::Traits::is_simd_vector<T>::value inline decltype(auto) fma(T &&x,
                                                                        T &&y,
                                                                        T &&z) {
    return Vc::fma(std::forward<T>(x), std::forward<T>(y), std::forward<T>(z));
}
/// @}

}  // namespace math
#elif defined(__CUDA_ARCH__)
namespace math {
using std::abs;

DETRAY_DEVICE inline float asin(float i) {
    return ::asinf(i);
}

DETRAY_DEVICE inline double asin(double i) {
    return ::asin(i);
}

DETRAY_DEVICE inline float atan(float i) {
    return ::atanf(i);
}

DETRAY_DEVICE inline double atan(double i) {
    return ::atan(i);
}

DETRAY_DEVICE inline float ceil(float i) {
    return ::ceilf(i);
}

DETRAY_DEVICE inline double ceil(double i) {
    return ::ceil(i);
}

DETRAY_DEVICE inline float copysign(float i, float j) {
    return ::copysignf(i, j);
}

DETRAY_DEVICE inline double copysign(double i, double j) {
    return ::copysign(i, j);
}

DETRAY_DEVICE inline float cos(float i) {
    return ::cosf(i);
}

DETRAY_DEVICE inline double cos(double i) {
    return ::cos(i);
}

DETRAY_DEVICE inline float exp(float i) {
    return ::expf(i);
}

DETRAY_DEVICE inline double exp(double i) {
    return ::exp(i);
}

DETRAY_DEVICE inline float fabs(float i) {
    return ::fabsf(i);
}

DETRAY_DEVICE inline double fabs(double i) {
    return ::fabs(i);
}

DETRAY_DEVICE inline float fma(float i, float j, float k) {
    return ::fmaf(i, j, k);
}

DETRAY_DEVICE inline double fma(double i, double j, double k) {
    return ::fma(i, j, k);
}

DETRAY_DEVICE inline float log(float i) {
    return ::logf(i);
}

DETRAY_DEVICE inline double log(double i) {
    return ::log(i);
}

DETRAY_DEVICE inline float log10(float i) {
    return ::log10f(i);
}

DETRAY_DEVICE inline double log10(double i) {
    return ::log10(i);
}

template <std::integral T>
DETRAY_DEVICE inline auto min(T i, T j) {
    return std::min(i, j);
}

DETRAY_DEVICE inline float min(float i, float j) {
    return ::min(i, j);
}

DETRAY_DEVICE inline double min(double i, double j) {
    return ::min(i, j);
}

template <std::integral T>
DETRAY_DEVICE inline auto max(T i, T j) {
    return std::max(i, j);
}

DETRAY_DEVICE inline float max(float i, float j) {
    return ::max(i, j);
}

DETRAY_DEVICE inline double max(double i, double j) {
    return ::max(i, j);
}

DETRAY_DEVICE inline float pow(float i, float p) {
    return ::powf(i, p);
}

DETRAY_DEVICE inline double pow(double i, double p) {
    return ::pow(i, p);
}

DETRAY_DEVICE inline auto signbit(float i) {
    return ::signbit(i);
}

DETRAY_DEVICE inline auto signbit(double i) {
    return ::signbit(i);
}

DETRAY_DEVICE inline float sin(float i) {
    return ::sinf(i);
}

DETRAY_DEVICE inline double sin(double i) {
    return ::sin(i);
}

DETRAY_DEVICE inline float sqrt(float i) {
    return ::sqrtf(i);
}

DETRAY_DEVICE inline double sqrt(double i) {
    return ::sqrt(i);
}

DETRAY_DEVICE inline float tan(float i) {
    return ::tanf(i);
}

DETRAY_DEVICE inline double tan(double i) {
    return ::tan(i);
}
}  // namespace math
#else
namespace math = std;
#endif  // SYCL

namespace detail {

using math::copysign;
using math::signbit;

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
    return math::copysign(mag, sgn);
#endif
}

/// Gets the signbit from a variable
template <typename scalar_t>
DETRAY_HOST_DEVICE inline bool signbit(scalar_t arg) {
#if defined(__CUDACC__)
    return signbit(arg);
#elif !defined(__CUDACC__)
    return math::signbit(arg);
#endif
}*/

}  // namespace detail

}  // namespace detray
