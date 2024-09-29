/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <cassert>
#include <concepts>
#include <iterator>
#include <type_traits>

/// Reimplement some usefull iterator functionality with __host__ __device__
/// quialifiers.
///
/// Adapted from libstdc++ e.g. @see
/// https://github.com/gcc-mirror/gcc/blob/16e2427f50c208dfe07d07f18009969502c25dc8/libstdc%2B%2B-v3/include/bits/stl_iterator_base_funcs.h
namespace detray::ranges {

// Reuse iterator tags for the host
using input_iterator_tag = std::input_iterator_tag;
using output_iterator_tag = std::output_iterator_tag;
using forward_iterator_tag = std::forward_iterator_tag;
using bidirectional_iterator_tag = std::bidirectional_iterator_tag;
using random_access_iterator_tag = std::random_access_iterator_tag;

// Iterator categories
template <class I>
inline constexpr bool input_iterator_v =
    std::is_base_of_v<detray::ranges::input_iterator_tag,
                      typename std::iterator_traits<I>::iterator_category>;

template <class I>
inline constexpr bool output_iterator_v =
    std::is_base_of_v<detray::ranges::output_iterator_tag,
                      typename std::iterator_traits<I>::iterator_category>;

template <class I>
inline constexpr bool forward_iterator_v =
    std::is_base_of_v<detray::ranges::forward_iterator_tag,
                      typename std::iterator_traits<I>::iterator_category>;

template <class I>
inline constexpr bool bidirectional_iterator_v =
    std::is_base_of_v<detray::ranges::bidirectional_iterator_tag,
                      typename std::iterator_traits<I>::iterator_category>;

template <class I>
inline constexpr bool random_access_iterator_v =
    std::is_base_of_v<detray::ranges::random_access_iterator_tag,
                      typename std::iterator_traits<I>::iterator_category>;

// Iterator concepts
template <class I>
concept input_iterator = detray::ranges::input_iterator_v<I>;

template <class I>
concept output_iterator = detray::ranges::output_iterator_v<I>;

template <class I>
concept forward_iterator = detray::ranges::forward_iterator_v<I>;

template <class I>
concept bidirectional_iterator = detray::ranges::bidirectional_iterator_v<I>;

template <class I>
concept random_access_iterator = detray::ranges::random_access_iterator_v<I>;

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
// TODO: ranges::ssize;
using std::data;
using std::empty;
// TODO: ranges::cdata
/// @}

///  @brief Reimplement std::distance.
/// @{
// used by every iterator up to and including bidirectional iterators
template <detray::ranges::input_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr
    typename std::iterator_traits<iterator_t>::difference_type
    distance_impl(iterator_t first, iterator_t last,
                  detray::ranges::input_iterator_tag) {
    typename std::iterator_traits<iterator_t>::difference_type d{0};
    // simply count
    while (first != last) {
        ++first;
        ++d;
    }
    return d;
}

// random access iterators specialization
template <detray::ranges::random_access_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr
    typename std::iterator_traits<iterator_t>::difference_type
    distance_impl(iterator_t first, iterator_t last,
                  detray::ranges::random_access_iterator_tag) {
    // use operator-
    return last - first;
}

template <detray::ranges::input_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr
    typename std::iterator_traits<iterator_t>::difference_type
    distance(iterator_t first, iterator_t last) {
    return distance_impl(
        first, last,
        typename std::iterator_traits<iterator_t>::iterator_category{});
}
/// @}

///  @brief Reimplement std::advance.
/// @{
template <detray::ranges::input_iterator iterator_t, typename dist_t>
DETRAY_HOST_DEVICE constexpr void advance_impl(
    iterator_t& itr, dist_t d, detray::ranges::input_iterator_tag) {
    static_assert(std::is_integral_v<dist_t>);
    assert(d > 0);
    // simply count
    while (d--) {
        ++itr;
    }
}

// random access iterators specialization
template <detray::ranges::bidirectional_iterator iterator_t, typename dist_t>
DETRAY_HOST_DEVICE constexpr void advance_impl(
    iterator_t& itr, dist_t d, detray::ranges::bidirectional_iterator_tag) {
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
template <detray::ranges::random_access_iterator iterator_t, typename dist_t>
DETRAY_HOST_DEVICE constexpr void advance_impl(
    iterator_t& itr, dist_t d, detray::ranges::random_access_iterator_tag) {
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

template <detray::ranges::input_iterator iterator_t, typename dist_t>
DETRAY_HOST_DEVICE constexpr void advance(iterator_t& itr, dist_t d) {
    return advance_impl(
        itr, d, typename std::iterator_traits<iterator_t>::iterator_category{});
}
/// @}

///  @brief Reimplement std::next and std::prev.
/// @{
template <detray::ranges::input_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr iterator_t next(
    iterator_t itr,
    typename std::iterator_traits<iterator_t>::difference_type d = 1) {
    detray::ranges::detail::advance(itr, d);
    return itr;
}

template <bidirectional_iterator iterator_t>
DETRAY_HOST_DEVICE constexpr iterator_t prev(
    iterator_t itr,
    typename std::iterator_traits<iterator_t>::difference_type d = 1) {
    detray::ranges::detail::advance(itr, -d);
    return itr;
}
/// @}

}  // namespace detail

}  // namespace detray::ranges
