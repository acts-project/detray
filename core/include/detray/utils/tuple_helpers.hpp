/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/type_traits.hpp"

// Thrust include(s)
#include <thrust/tuple.h>

// System include(s)
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

namespace detray::detail {

/// get function accessor for either std::tuple or thrust::tuple
///
/// usage example:
/// detail::get<0>(tuple)
/// @{
using std::get;
using thrust::get;

/// Retrieve an element from a thrust tuple by value. No perfect forwarding for
/// composite types like tuple_t<value_types...>
template <typename query_t, template <typename...> class tuple_t,
          class... value_types,
          std::enable_if_t<std::is_same_v<tuple_t<value_types...>,
                                          thrust::tuple<value_types...>>,
                           bool> = true>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    const tuple_t<value_types...>& tuple) noexcept {
    return thrust::get<get_type_pos_v<query_t, value_types...>>(
        std::forward<const tuple_t<value_types...>>(tuple));
}

/// Retrieve an element from a thrust tuple by value. No perfect forwarding for
/// composite types like tuple_t<value_types...>
template <typename query_t, template <typename...> class tuple_t,
          class... value_types,
          std::enable_if_t<std::is_same_v<tuple_t<value_types...>,
                                          thrust::tuple<value_types...>>,
                           bool> = true>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    tuple_t<value_types...>& tuple) noexcept {
    return thrust::get<get_type_pos_v<query_t, value_types...>>(
        std::forward<tuple_t<value_types...>>(tuple));
}
/// @}

/// tuple_element for either std::tuple or thrust::tuple
///
/// usage example:
/// detail::tuple_element< int, tuple_t >::type
/// @{
template <int N, class T, typename Enable = void>
struct tuple_element;

// std::tuple
template <int N, template <typename...> class tuple_t, class... value_types>
struct tuple_element<N, tuple_t<value_types...>,
                     typename std::enable_if_t<
                         std::is_same_v<tuple_t<value_types...>,
                                        std::tuple<value_types...>> == true>>
    : std::tuple_element<N, tuple_t<value_types...>> {};

// thrust::tuple
template <int N, template <typename...> class tuple_t, class... value_types>
struct tuple_element<N, tuple_t<value_types...>,
                     typename std::enable_if_t<
                         std::is_same_v<tuple_t<value_types...>,
                                        thrust::tuple<value_types...>> == true>>
    : thrust::tuple_element<N, tuple_t<value_types...>> {};
/// @}

/// tuple_size for either std::tuple or thrust::tuple
///
/// usage example:
/// detail::tuple_size< tuple_t >::value
/// @{
template <class T, typename Enable = void>
struct tuple_size;

// std::tuple
template <template <typename...> class tuple_t, class... value_types>
struct tuple_size<tuple_t<value_types...>,
                  typename std::enable_if_t<
                      std::is_same_v<tuple_t<value_types...>,
                                     std::tuple<value_types...>> == true>>
    : std::tuple_size<tuple_t<value_types...>> {};

// thrust::tuple
template <template <typename...> class tuple_t, class... value_types>
struct tuple_size<tuple_t<value_types...>,
                  typename std::enable_if_t<
                      std::is_same_v<tuple_t<value_types...>,
                                     thrust::tuple<value_types...>> == true>>
    : thrust::tuple_size<tuple_t<value_types...>> {};
/// @}

/// make_tuple for either std::tuple or thrust::tuple
/// users have to specifiy tuple_t for detail::make_tuple
///
/// usage example
/// detail::make_tuple<tuple_t>(args...)
/// @{
template <class T>
struct unwrap_refwrapper {
    using type = T;
};

template <class T>
struct unwrap_refwrapper<std::reference_wrapper<T>> {
    using type = T&;
};

template <class T>
using unwrap_decay_t =
    typename unwrap_refwrapper<typename std::decay<T>::type>::type;

// make_tuple for std::tuple
template <template <typename...> class tuple_t, class... value_types,
          std::enable_if_t<std::is_same_v<tuple_t<value_types...>,
                                          std::tuple<value_types...>>,
                           bool> = true>
DETRAY_HOST_DEVICE inline constexpr std::tuple<unwrap_decay_t<value_types>...>
make_tuple(value_types&&... args) {
    return std::make_tuple(args...);
}

// make_tuple for thrust::tuple
template <template <typename...> class tuple_t, class... value_types,
          std::enable_if_t<std::is_same_v<tuple_t<value_types...>,
                                          thrust::tuple<value_types...>>,
                           bool> = true>
DETRAY_HOST_DEVICE inline constexpr thrust::tuple<
    unwrap_decay_t<value_types>...>
make_tuple(value_types&&... args) {
    return thrust::make_tuple(args...);
}
/// @}

}  // namespace detray::detail