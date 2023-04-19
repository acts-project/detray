/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/sort.hpp"

// Thrust include(s)
#if defined(__CUDACC__)
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>
#include <thrust/find.h>
#include <thrust/sort.h>
#endif

// System include(s)
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <utility>

namespace detray::detail {

/// Composes a floating point value with the magnitude of @param mag and the
/// sign of @param sgn
template <typename scalar_t>
DETRAY_HOST_DEVICE inline scalar_t copysign(scalar_t mag, scalar_t sgn) {
#if defined(__CUDACC__)
    if constexpr (std::is_same_v<scalar_t, float>) {
        return copysignf(mag, sgn);
    } else {
        return copysign(mag, sgn);
    }
#elif !defined(__CUDACC__)
    return std::copysign(mag, sgn);
#endif
}

/// @brief sequential (single thread) sort function
template <class RandomIt>
DETRAY_HOST_DEVICE inline void sequential_sort(RandomIt first, RandomIt last) {
#if defined(__CUDACC__) || defined(CL_SYCL_LANGUAGE_VERSION) || \
    defined(SYCL_LANGUAGE_VERSION)
    selection_sort(first, last);
#else
    std::sort(first, last);
#endif
}

/// @brief sequential sort for host/devcie (single thread) - custom comparison
template <class RandomIt, class Compare>
DETRAY_HOST_DEVICE inline void sequential_sort(RandomIt first, RandomIt last,
                                               Compare&& comp) {
#if defined(__CUDACC__)
    thrust::sort(thrust::seq, first, last, comp);
#elif !defined(__CUDACC__)
    std::sort(first, last, std::forward<Compare>(comp));
#endif
}

/// @brief find_if implementation for host/devcie (single thread)
template <class RandomIt, class Predicate>
DETRAY_HOST_DEVICE inline auto find_if(RandomIt first, RandomIt last,
                                       Predicate&& comp) {
    for (RandomIt i = first; i != last; ++i) {
        if (comp(*i)) {
            return i;
        }
    }

    return last;
}

/// @brief lower_bound implementation for device
template <class ForwardIt, typename Value>
DETRAY_HOST_DEVICE inline auto lower_bound(ForwardIt first, ForwardIt last,
                                           Value value) {
#if defined(__CUDACC__)
    return thrust::lower_bound(thrust::seq, first, last, value);
#elif !defined(__CUDACC__)
    return std::lower_bound(first, last, value);
#endif
}

/// @brief upper_bound implementation for device
template <class ForwardIt, typename Value>
DETRAY_HOST_DEVICE inline auto upper_bound(ForwardIt first, ForwardIt last,
                                           Value value) {
#if defined(__CUDACC__)
    return thrust::upper_bound(thrust::seq, first, last, value);
#elif !defined(__CUDACC__)
    return std::upper_bound(first, last, value);
#endif
}

}  // namespace detray::detail
