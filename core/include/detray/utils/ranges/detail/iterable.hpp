/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <cstddef>
#include <iterator>
#include <type_traits>

namespace detray::ranges::detail {

/// @brief Detail interface of an iterable type.
///
/// Provides @c begin() and @c end() functions of different container
/// implementations. These can be used as a range.
///
/// Base case: The type is not iterable
template <typename, typename = void>
struct iterable : public std::false_type {};

/// @brief Iterable specialization.
///
/// The type is iterable by simply calling its begin and end member functions.
/// In case of @c vecmem::device_vector the iterator is a pointer type.
// TODO: Use concepts in the future. Right now, the existance of 'begin()' is
// used to infer on the existence of 'end()' etc.
template <typename T>
struct iterable<
    T, std::enable_if_t<std::is_class_v<std::decay_t<decltype(std::begin(
                            std::declval<T&>()))>> or
                            std::is_pointer_v<std::decay_t<decltype(std::begin(
                                std::declval<T&>()))>>,
                        void>> : public std::true_type {

    using iterator_type =
        std::decay_t<decltype(std::begin(std::declval<T&>()))>;
    using const_iterator_type =
        std::decay_t<decltype(std::cbegin(std::declval<T&>()))>;

    // non-const

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline constexpr auto begin(T& iterable) noexcept -> iterator_type {
        return std::begin(iterable);
    }

    /// @returns sentinel, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline constexpr auto end(T& iterable) noexcept -> iterator_type {
        return std::end(iterable);
    }

    // const

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline constexpr auto begin(const T& iterable) noexcept
        -> const_iterator_type {
        return std::cbegin(iterable);
    }

    /// @returns sentinel, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline constexpr auto end(const T& iterable) noexcept
        -> const_iterator_type {
        return std::cend(iterable);
    }
};

/// @brief Iterable specialization for single object outside of a container.
///
/// Provides a @c begin() and @c end() implementation that can be used with
/// iterator aware code.
// TODO: Use concepts in the future. Right now, the existance of 'begin()' is
// used to infer on the compliance with the rest of the iterator interface
template <typename T>
struct iterable<
    T, std::enable_if_t<not std::is_class_v<typename T::value_type>, bool>>
    : public std::true_type {

    using iterator_type = T*;
    using const_iterator_type = const T*;

    // non-const

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline constexpr auto begin(T& iterable) noexcept -> iterator_type {
        return &iterable;
    }

    /// @returns sentinel, which is the position behind the wrapped value.
    /// @note protected against being modified.
    DETRAY_HOST_DEVICE
    static inline constexpr auto end(T& iterable) noexcept -> iterator_type {
        return &iterable + std::ptrdiff_t{1};
    }

    // const

    /// @returns start position, which is at the wrapped value.
    DETRAY_HOST_DEVICE
    static inline constexpr auto begin(const T& iterable) noexcept
        -> const_iterator_type {
        return &iterable;
    }

    /// @returns sentinel, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline constexpr auto end(const T& iterable) noexcept
        -> const_iterator_type {
        return &iterable + std::ptrdiff_t{1};
    }
};

}  // namespace detray::ranges::detail