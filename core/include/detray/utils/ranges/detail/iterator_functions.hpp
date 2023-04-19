/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Thrust include(s)
#if defined(__CUDACC__)
#include <thrust/iterator/iterator_categories.h>
#endif

// Project include(s)
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <cassert>
#include <iterator>
#include <type_traits>

/// Reimplement some usefull iterator functionality with __host__ __device__
/// quialifiers.
///
/// Adapted from libstdc++ e.g. @see
/// https://github.com/gcc-mirror/gcc/blob/16e2427f50c208dfe07d07f18009969502c25dc8/libstdc%2B%2B-v3/include/bits/stl_iterator_base_funcs.h
namespace detray::ranges {

// Define iterator tags for host and device
/*#if defined(__CUDACC__)
using input_iterator_tag = thrust::input_device_iterator_tag;
using output_iterator_tag = thrust::output_device_iterator_tag;
using forward_iterator_tag = thrust::forward_device_iterator_tag;
using bidirectional_iterator_tag = thrust::bidirectional_device_iterator_tag;
using random_access_iterator_tag = thrust::random_access_device_iterator_tag;
#elif !defined(__CUDACC__)*/
// reuse iterator tags for the host
using input_iterator_tag = std::input_iterator_tag;
using output_iterator_tag = std::output_iterator_tag;
using forward_iterator_tag = std::forward_iterator_tag;
using bidirectional_iterator_tag = std::bidirectional_iterator_tag;
using random_access_iterator_tag = std::random_access_iterator_tag;
// #endif

namespace detail {
/// Simply import the std versions of basic iterator functionality where
/// possible
/// @{
using std::begin;
using std::cbegin;
using std::cend;
using std::crbegin;
using std::crend;
using std::end;
using std::rbegin;
using std::rend;

using std::size;
// ranges::ssize;
using std::data;
using std::empty;
// ranges::cdata
/// @}

///  @brief Reimplement std::distance.
/// @{
// used by every iterator up to and including bidirectional iterators
template <class input_iterator_t>
DETRAY_HOST_DEVICE inline constexpr
    typename std::iterator_traits<input_iterator_t>::difference_type
    distance_impl(input_iterator_t first, input_iterator_t last,
                  detray::ranges::input_iterator_tag) {
    typename std::iterator_traits<input_iterator_t>::difference_type d{0};
    // simply count
    while (first != last) {
        ++first;
        ++d;
    }
    return d;
}

// random access iterators specialization
template <class random_access_iterator_t>
DETRAY_HOST_DEVICE inline constexpr
    typename std::iterator_traits<random_access_iterator_t>::difference_type
    distance_impl(random_access_iterator_t first, random_access_iterator_t last,
                  detray::ranges::random_access_iterator_tag) {
    // use operator-
    return last - first;
}

template <
    class iterator_t,
    std::enable_if_t<std::is_base_of_v<detray::ranges::input_iterator_tag,
                                       typename std::iterator_traits<
                                           iterator_t>::iterator_category>,
                     bool> = true>
DETRAY_HOST_DEVICE inline constexpr
    typename std::iterator_traits<iterator_t>::difference_type
    distance(iterator_t first, iterator_t last) {
    return distance_impl(
        first, last,
        typename std::iterator_traits<iterator_t>::iterator_category{});
}
/// @}

///  @brief Reimplement std::advance.
/// @{
template <class input_iterator_t, typename dist_t>
DETRAY_HOST_DEVICE inline constexpr void advance_impl(
    input_iterator_t& itr, dist_t d, detray::ranges::input_iterator_tag) {
    static_assert(std::is_integral_v<dist_t>);
    assert(d > 0);
    // simply count
    while (d--) {
        ++itr;
    }
}

// random access iterators specialization
template <class bidirectional_iterator_t, typename dist_t>
DETRAY_HOST_DEVICE inline constexpr void advance_impl(
    bidirectional_iterator_t& itr, dist_t d,
    detray::ranges::bidirectional_iterator_tag) {
    static_assert(std::is_integral_v<dist_t>);
    if (d > 0) {
        while (d--) {
            ++itr;
        }
    } else {
        while (d++) {
            --itr;
        }
    }
}

// random access iterators specialization
template <class random_access_iterator_t, typename dist_t>
DETRAY_HOST_DEVICE inline constexpr void advance_impl(
    random_access_iterator_t& itr, dist_t d,
    detray::ranges::random_access_iterator_tag) {
    static_assert(std::is_integral_v<dist_t>);
    if (d == static_cast<dist_t>(1)) {
        ++itr;
    } else {
        if constexpr (std::is_signed_v<dist_t>) {
            if (d == static_cast<dist_t>(-1)) {
                --itr;
                return;
            }
        }
        itr += d;
    }
}

template <
    class iterator_t, typename dist_t,
    std::enable_if_t<std::is_base_of_v<detray::ranges::input_iterator_tag,
                                       typename std::iterator_traits<
                                           iterator_t>::iterator_category>,
                     bool> = true>
DETRAY_HOST_DEVICE inline constexpr void advance(iterator_t& itr, dist_t d) {
    return advance_impl(
        itr, d, typename std::iterator_traits<iterator_t>::iterator_category{});
}
/// @}

///  @brief Reimplement std::next and std::prev.
/// @{
template <class input_iterator_t>
DETRAY_HOST_DEVICE inline constexpr input_iterator_t next(
    input_iterator_t itr,
    typename std::iterator_traits<input_iterator_t>::difference_type d = 1) {
    detray::ranges::detail::advance(itr, d);
    return itr;
}

template <
    class bidirectional_iterator_t,
    std::enable_if_t<
        std::is_base_of_v<detray::ranges::bidirectional_iterator_tag,
                          typename std::iterator_traits<
                              bidirectional_iterator_t>::iterator_category>,
        bool> = true>
DETRAY_HOST_DEVICE inline constexpr bidirectional_iterator_t prev(
    bidirectional_iterator_t itr,
    typename std::iterator_traits<bidirectional_iterator_t>::difference_type d =
        1) {
    detray::ranges::detail::advance(itr, -d);
    return itr;
}
/// @}

}  // namespace detail

}  // namespace detray::ranges
