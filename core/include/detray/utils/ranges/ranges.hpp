/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <type_traits>

#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/detail/iterable.hpp"
#include "detray/utils/ranges/detail/iterator_traits.hpp"

namespace detray::ranges {

/// Shorthand to test a type for being iterable.
template <typename T>
inline constexpr bool is_iterable_v = detail::iterable<T>::value;

/// @brief Detail interface of an iterable type.
///
/// Base case: The type is not iterable
template <typename T, typename = void>
struct range : public std::false_type {};

/// @brief Check if a type is iterable and extract iterator types.
///
/// A range has a begin() and an end() member function, thus becoming iterable.
/// It wraps some kind of iterable: this can be a container which is either
/// owned by the range or not, or it can be another instace of @c iterable, e.g.
/// another range.
template <typename T>
struct range<T, std::enable_if_t<is_iterable_v<T>, void>>
    : public std::true_type {

    /// Extract begin() and end() functions from [detray] containers
    using iterable_type = detail::iterable<T>;

    /// Compliance with std::ranges typedefs
    using iterator_t = typename iterable_type::iterator_type;
    using const_iterator_t = typename iterable_type::const_iterator_type;
    using range_value_t =
        typename detail::iterator_traits<iterator_t>::value_type;
    using range_size_t = std::size_t;
    using range_difference_t =
        typename detail::iterator_traits<iterator_t>::difference_type;
    using range_reference_t =
        typename detail::iterator_traits<iterator_t>::reference;
    using range_const_reference_t = const range_reference_t;
    using range_rvalue_reference_t = range_value_t&&;
};

/// Check if a type @c T is a range and get the @c begin function
template <class T>
auto begin(T&& iterable) {
    return range<T>::iterable_type::begin(iterable);
}

/// Check if a type @c T is a range and get the @c end function
template <class T>
auto end(T&& iterable) {
    return range<T>::iterable_type::end(iterable);
}

/// Check if a type @c T is a range and get the @c rbegin function
template <class T>
auto rbegin(T&& iterable) {
    return range<T>::iterable_type::rbegin(iterable);
}

/// Check if a type @c T is a range and get the @c rend function
template <class T>
auto rend(T&& iterable) {
    return range<T>::iterable_type::rend(iterable);
}

/// Check if a type @c T is a range and get the @c cbegin function
template <class T>
auto cbegin(T&& iterable) {
    return range<T>::iterable_type::cbegin(iterable);
}

/// Check if a type @c T is a range and get the @c cend function
template <class T>
auto cend(T&& iterable) {
    return range<T>::iterable_type::cend(iterable);
}

/// Check if a type @c T is a range and get the @c crbegin function
template <class T>
auto crbegin(T&& iterable) {
    return range<T>::iterable_type::crbegin(iterable);
}

/// Check if a type @c T is a range and get the @c crend function
template <class T>
auto crend(T&& iterable) {
    return range<T>::iterable_type::crend(iterable);
}

}  // namespace detray::ranges

// https://en.cppreference.com/w/cpp/ranges/view
namespace detray::views {

/// Tags a type as a view
struct base_view {};

/// Defines a detray 'view'
template <typename view_impl_t>
struct view_interface {

    view_interface() = default;

    DETRAY_HOST_DEVICE
    auto empty() const { return (size() == 0); }

    DETRAY_HOST_DEVICE
    auto data() const { return m_impl_ptr->begin(); }

    DETRAY_HOST_DEVICE
    auto size() const {
        return std::distance(m_impl_ptr->end(), m_impl_ptr->begin());
    }

    DETRAY_HOST_DEVICE
    auto front() { return *(m_impl_ptr->begin()); }

    DETRAY_HOST_DEVICE
    auto back() { return *(m_impl_ptr->end()); }

    /// @return element at position i, relative to iterator range.
    DETRAY_HOST_DEVICE
    inline auto& operator[](const dindex i) {
        return *(m_impl_ptr->begin() + i);
    }

    /// @return element at position i, relative to iterator range - const
    DETRAY_HOST_DEVICE
    inline const auto& operator[](const dindex i) const {
        return *(m_impl_ptr->begin() + i);
    }

    const view_impl_t* m_impl_ptr{static_cast<const view_impl_t*>(this)};
};

/*template<template<typename> class T, typename iterable_t>
inline constexpr bool enable_view
    = std::is_base_of_v<base_view, T<iterable_t>> or
        std::is_base_of_v<view_interface<T<iterable_t>>, T<iterable_t>>;

// https://en.cppreference.com/w/cpp/ranges/owning_view
template<template range_t>
struct owning_view<range_t> : public view_interface<owning_view<range_t>> {
    owning_view() = delete;

    owning_view(range_t && data) : m_data(std::move(data)) {}

    /// @returns start position, which is at the wrapped value.
    template<typename... Args>
    DETRAY_HOST_DEVICE
     inline auto begin() const -> iterator_type {
        return ranges::begin(m_data);
    }

    /// @returns end position, which is the position behind the wrapped value.
    DETRAY_HOST_DEVICE
    static inline auto end() const -> iterator_type {
        return ranges::end(m_data);
    }

    range_t m_data;
};

template<template range_t>
struct non_owning_view<range_t>
    : public view_interface<non_owning_view<range_t>> {

    non_owning_view() = default;

    non_owning_view(const range_t & data) : m_data(data) {}

    const range_t & m_data;
};*/

}  // namespace detray::views
