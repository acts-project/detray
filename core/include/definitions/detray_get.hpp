/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <array>
#include <tuple>
#include <utility>

#if defined(__CUDACC__)
#include <thrust/tuple.h>
#endif

#include "definitions/cuda_qualifiers.hpp"

namespace detray {

namespace detail {

/** This header defines how to get elements from detray container and range
 * types in generic parts of the code.
 */

#if defined(__CUDACC__)
template <unsigned int index, typename... value_types>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    const thrust::tuple<value_types...> &tuple) {
    return thrust::get<index>(tuple);
}

template <unsigned int index, typename... value_types>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    thrust::tuple<value_types...> &tuple) {
    return thrust::get<index>(tuple);
}

template <unsigned int index, typename... value_types>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    const thrust::pair<value_types...> &pair) {
    return thrust::get<index>(pair);
}

template <unsigned int index, typename... value_types>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    thrust::pair<value_types...> &pair) {
    return thrust::get<index>(pair);
}
#endif

template <unsigned int index, typename... value_types>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    const std::tuple<value_types...> &tuple) {
    return std::get<index>(tuple);
}

template <unsigned int index, typename... value_types>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    std::tuple<value_types...> &tuple) {
    return std::get<index>(tuple);
}

/** In case the index is not known at compile time. */
template <unsigned int current_index = 0, typename... value_types>
DETRAY_HOST_DEVICE inline decltype(auto) get(
    const std::tuple<value_types...> &tuple, unsigned int index) {
    // Check index against constexpr
    if (current_index == index) {
        return std::get<current_index>(tuple);
    }
    // Next tuple index
    if constexpr (current_index < std::tuple_size_v<tuple> - 1) {
        return get<current_index + 1>(tuple, index);
    }
}

template <unsigned int current_index = 0, typename... value_types>
DETRAY_HOST_DEVICE inline decltype(auto) get(std::tuple<value_types...> &tuple,
                                             unsigned int index) {
    // Check index against constexpr
    if (current_index == index) {
        return std::get<current_index>(tuple);
    }
    // Next tuple index
    if constexpr (current_index < std::tuple_size_v<tuple> - 1) {
        return get<current_index + 1>(tuple, index);
    }
}

template <unsigned int index, typename... value_types>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    const std::pair<value_types...> &pair) {
    return std::get<index>(pair);
}

template <unsigned int index, typename... value_types>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    std::pair<value_types...> &pair) {
    return std::get<index>(pair);
}

template <unsigned int index, typename value_type, unsigned int N>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    std::array<value_type, N> &array) {
    return std::get<index>(array);
}

template <unsigned int index, typename value_type, unsigned int N>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    const std::array<value_type, N> &array) {
    return std::get<index>(array);
}

/** Put new link types here */

}  // namespace detail
}  // namespace detray