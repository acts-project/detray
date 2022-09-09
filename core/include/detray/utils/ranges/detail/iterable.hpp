/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <cstddef>
#include <iterator>
#include <type_traits>

#include "detray/definitions/qualifiers.hpp"

namespace detray::ranges::detail {

/// @brief Detail interface of an iterable type.
///
/// Provides @c begin() and @c end() functions to different container
/// implementations. These can be used to define a range.
///
/// Base case: The type is not iterable
template <typename, typename = void>
struct iterable : public std::false_type {};

/// @brief Iterable specialization.
///
/// The type is iterable by simply calling its begin and end member functions.
// TODO: Use concepts in the future. Right now, the existance of 'begin()' is
// used to infer on the compliance with the
template <typename T>
struct iterable<
    T, std::enable_if_t<
           std::is_class_v<decltype(std::begin(std::declval<T>()))>, void>>
    : public std::true_type {
    using iterator_type = decltype(std::begin(std::declval<T>()));
    using reverse_iterator_type = decltype(std::rbegin(std::declval<T>()));
    using const_iterator_type = decltype(std::cbegin(std::declval<T>()));
    using const_reverse_iterator_type =
        decltype(std::crbegin(std::declval<T>()));

    // non-const

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto begin(T& iterable) -> iterator_type {
        return iterable.begin();
    }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto end(T& iterable) -> iterator_type {
        return iterable.end();
    }

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto rbegin(T& iterable) -> reverse_iterator_type {
        return iterable.rbegin();
    }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto rend(T& iterable) -> reverse_iterator_type {
        return iterable.rend();
    }

    // const

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto cbegin(const T& iterable) -> const_iterator_type {
        return iterable.cbegin();
    }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto cend(const T& iterable) -> const_iterator_type {
        return iterable.cend();
    }

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto crbegin(const T& iterable)
        -> const_reverse_iterator_type {
        return iterable.crbegin();
    }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto crend(const T& iterable) -> const_reverse_iterator_type {
        return iterable.crend();
    }
};

/// @brief Iterable specialization.
///
/// The type is iterable by simply calling its begin and end member functions.
// TODO: Use concepts in the future. Right now, the existance of 'begin()' is
// used to infer on the compliance with the rest of the iterator interface
/*template<typename T>
struct iterable<T, std::enable_if_t<std::is_pointer_v<T>, bool>>
    : public std::true_type {
    using iterator_type = T;
    using reverse_iterator_type = T;
    using const_iterator_type = const T;
    using const_reverse_iterator_type = const T;

    // non-const

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto begin(T & iterable) -> iterator_type {
        return iterable;
    }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto end(T & iterable, const std::ptrdiff_t n) ->
iterator_type { return iterable + n;
    }

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto rbegin(T & iterable, const std::ptrdiff_t n)
        -> reverse_iterator_type {
        return iterable + n;
    }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto rend(T & iterable)
        -> reverse_iterator_type {
        return iterable;
    }

    // const

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto cbegin(const T & iterable) -> const_iterator_type {
        return iterable;
    }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto cend(const T & iterable, const std::ptrdiff_t n) ->
const_iterator_type { return iterable + n;
    }

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto crbegin(const T & iterable)
        -> const_reverse_iterator_type {
        return iterable;
    }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto crend(const T & iterable, const std::ptrdiff_t n)
        -> const_reverse_iterator_type {
        return iterable + n;
    }
};*/

/// @brief Iterable specialization.
///
/// Specialization for detray specific cases, e.g. the tuple container's groups
/*template<typename T, typename itr_t = ...>
struct iterable : public std::true_type { };
*/

}  // namespace detray::ranges::detail