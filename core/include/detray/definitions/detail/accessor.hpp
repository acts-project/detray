/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <thrust/tuple.h>

#include <array>
#include <tuple>
#include <utility>

#include "detray/definitions/detray_qualifiers.hpp"

namespace detray {

/** define tuple type namespace for host (std) and device (thrust) compiler
 * **/
#if defined(__CUDACC__)
namespace vtuple = thrust;
#else
namespace vtuple = std;
#endif

namespace detail {

/** tuple size
 */
template <class T>
struct tuple_size;

template <template <typename...> class tuple_type, class... value_types>
struct tuple_size<tuple_type<value_types...>>
    : std::integral_constant<std::size_t, sizeof...(value_types)> {};

/** thrust tuple accessor
 */

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

/** thrust pair accessor
 */

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

/** std tuple accessor
 */

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
    if constexpr (current_index <
                  std::tuple_size_v<std::tuple<value_types...>> - 1) {
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
    if constexpr (current_index <
                  std::tuple_size_v<std::tuple<value_types...>> - 1) {
        return get<current_index + 1>(tuple, index);
    }
}

/** std pair accessor
 */

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

/** std array accessor
 */

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