/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#if defined(__CUDACC__)
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>
#include <thrust/find.h>
#include <thrust/sort.h>
#endif

#include <thrust/tuple.h>

#include <algorithm>
#include <array>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

#include "detray/definitions/qualifiers.hpp"

namespace detray {

namespace detail {

/// Extract the first type from a parameter pack.
template <typename first_t = void, typename... other_t>
struct first_type {
    using type = first_t;
};

template <typename... TYPES>
using first_t = typename first_type<TYPES...>::type;

/// Extract the first value from an index pack.
template <std::size_t first_t, std::size_t... other_t>
struct first_idx {
    static constexpr std::size_t value = first_t;
};

template <std::size_t... TYPES>
inline constexpr std::size_t first_idx_v = first_idx<TYPES...>::value;

namespace {
struct empty_type {};

/// Unroll a parameter pack without using a tuple.
///
/// Returns the position of the type counted from the back!
template <typename query_t, typename first_t = empty_type,
          typename... remaining_types>
DETRAY_HOST_DEVICE constexpr std::size_t unroll_values() {
    if constexpr (not std::is_same_v<first_t, empty_type> and
                  not std::is_same_v<query_t, first_t>) {
        return unroll_values<query_t, remaining_types...>();
    }
    if constexpr (std::is_same_v<query_t, first_t>) {
        return sizeof...(remaining_types) + 1;
    }
    return std::numeric_limits<std::size_t>::max();
}

}  // anonymous namespace

/** get function accessor
 *
 *  usage example:
 *  detail::get<0>(tuple)
 */
using std::get;
using thrust::get;

template <std::size_t id, typename mask_store_t,
          std::enable_if_t<std::is_class_v<typename std::remove_reference_t<
                               mask_store_t>::container_type>,
                           bool> = true>
DETRAY_HOST_DEVICE constexpr auto get(mask_store_t&& mask_store) noexcept
    -> decltype(get<id>(std::forward<mask_store_t>(mask_store).get())) {
    return get<id>(std::forward<mask_store_t>(mask_store).get());
}

/// Retrieve an element from a thrust tuple by value. No perfect forwarding for
/// composite types like tuple_t<value_types...>
template <typename query_t, template <typename...> class tuple_t,
          class... value_types,
          std::enable_if_t<std::is_same_v<tuple_t<value_types...>,
                                          thrust::tuple<value_types...>>,
                           bool> = true>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const tuple_t<value_types...>& tuple) noexcept {
    return thrust::get<sizeof...(value_types) -
                       unroll_values<query_t, value_types...>()>(
        std::forward<const tuple_t<value_types...>>(tuple));
}

/// Retrieve an element from a thrust tuple by value. No perfect forwarding for
/// composite types like tuple_t<value_types...>
template <typename query_t, template <typename...> class tuple_t,
          class... value_types,
          std::enable_if_t<std::is_same_v<tuple_t<value_types...>,
                                          thrust::tuple<value_types...>>,
                           bool> = true>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    tuple_t<value_types...>& tuple) noexcept {
    return thrust::get<sizeof...(value_types) -
                       unroll_values<query_t, value_types...>()>(
        std::forward<tuple_t<value_types...>>(tuple));
}

/**
 *  tuple_element accessor
 *
 *  usage example:
 *  detail::tuple_element< int, tuple_t >::type
 */
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

/**
 *  tuple size accessor
 *
 *  usage example:
 *  detail::tuple_size< tuple_t >::value
 */
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

/**
 *  make_tuple accessor
 *  users have to specifiy tuple_t for detail::make_tuple
 *
 *  usage example
 *  detail::make_tuple<tuple_t>(args...)
 */

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
DETRAY_HOST_DEVICE constexpr std::tuple<unwrap_decay_t<value_types>...>
make_tuple(value_types&&... args) {
    return std::make_tuple(args...);
}

// make_tuple for thrust::tuple
template <template <typename...> class tuple_t, class... value_types,
          std::enable_if_t<std::is_same_v<tuple_t<value_types...>,
                                          thrust::tuple<value_types...>>,
                           bool> = true>
DETRAY_HOST_DEVICE constexpr thrust::tuple<unwrap_decay_t<value_types>...>
make_tuple(value_types&&... args) {
    return thrust::make_tuple(args...);
}

/**
 *  Trait class to figure out if a given type has a @c reserve(...) function
 */
template <typename T>
struct has_reserve {

    private:
    /// Function returning @c std::true_type for types that do have a @c
    /// reserve(...) function
    template <typename C>
    static constexpr auto check(C*) ->
        typename std::is_void<decltype(std::declval<C>().reserve(
            std::declval<typename C::size_type>()))>::type;
    /// Function returning @c std::false_type for types that fair the previous
    /// function
    template <typename>
    static constexpr std::false_type check(...);

    /// Declare the value type of this trait class
    typedef decltype(check<T>(nullptr)) type;

    public:
    /// Value of the check
    static constexpr bool value = type::value;
};

/// @name Functions calling or not calling reserve(...) based on whether it's
/// available
/// @{

template <typename T, std::enable_if_t<has_reserve<T>::value, bool> = true>
DETRAY_HOST_DEVICE void call_reserve(T& obj, std::size_t newsize) {
    obj.reserve(newsize);
}

template <typename T, std::enable_if_t<!(has_reserve<T>::value), bool> = true>
DETRAY_HOST_DEVICE void call_reserve(T& /*obj*/, std::size_t /*newsize*/) {}

/**
 *  sequential (single thread) sort function
 */

template <class RandomIt>
DETRAY_HOST_DEVICE void sequential_sort(RandomIt first, RandomIt last) {

#if defined(__CUDACC__)
    thrust::sort(thrust::seq, first, last);
#elif !defined(__CUDACC__)
    std::sort(first, last);
#endif
}

template <class RandomIt, class Compare>
DETRAY_HOST_DEVICE void sequential_sort(RandomIt first, RandomIt last,
                                        Compare comp) {
#if defined(__CUDACC__)
    thrust::sort(thrust::seq, first, last, comp);
#elif !defined(__CUDACC__)
    std::sort(first, last, comp);
#endif
}

/**
 *  find_if implementation for device
 */
template <class RandomIt, class Predicate>
DETRAY_HOST_DEVICE auto find_if(RandomIt first, RandomIt last, Predicate comp) {
#if defined(__CUDACC__)
    return thrust::find_if(thrust::seq, first, last, comp);
#elif !defined(__CUDACC__)
    return std::find_if(first, last, comp);
#endif
}

/// @brief lower_bound implementation for device
template <class ForwardIt, typename Value>
DETRAY_HOST_DEVICE auto lower_bound(ForwardIt first, ForwardIt last,
                                    Value value) {
#if defined(__CUDACC__)
    return thrust::lower_bound(thrust::seq, first, last, value);
#elif !defined(__CUDACC__)
    return std::lower_bound(first, last, value);
#endif
}

/// @brief upper_bound implementation for device
template <class ForwardIt, typename Value>
DETRAY_HOST_DEVICE auto upper_bound(ForwardIt first, ForwardIt last,
                                    Value value) {
#if defined(__CUDACC__)
    return thrust::upper_bound(thrust::seq, first, last, value);
#elif !defined(__CUDACC__)
    return std::upper_bound(first, last, value);
#endif
}

}  // namespace detail
}  // namespace detray
