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
#elif(IS_SOA)
#include <Vc/Vc>
#endif

// System include(s).
#include <cmath>

namespace detray {

/// Namespace to pick up math functions from
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
namespace math_ns = cl::sycl;
#else
namespace math_ns = std;
#endif  // SYCL

namespace detail {

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
