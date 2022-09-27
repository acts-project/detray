/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <iterator>
#include <type_traits>

namespace detray::ranges {

/// @brief Detail interface of an iterable type.
///
/// Base case: The type is not iterable
template <typename T, typename = void>
struct range : public std::false_type {
    using type = void;
};

/// @brief Check if a type is iterable and extract iterator type.
///
/// A range has a begin() and an end() member function, which is guranteed by
/// simply calling std::begin and std::end.
/// In case of @c vecmem::device_vector the iterator is a pointer type.
template <typename T>
struct range<
    T, std::enable_if_t<std::is_class_v<std::decay_t<decltype(std::begin(
                            std::declval<T&>()))>> or
                            std::is_pointer_v<std::decay_t<decltype(std::begin(
                                std::declval<T&>()))>>,
                        void>> : public std::true_type {
    using iterator_type =
        std::decay_t<decltype(std::begin(std::declval<T&>()))>;
    using const_iterator_type =
        std::decay_t<decltype(std::cbegin(std::declval<T&>()))>;
};

class base_view;

template <typename range_t>
class view_interface;

template <typename T>
inline constexpr bool enable_view =
    std::is_base_of_v<base_view, T> or std::is_base_of_v<view_interface<T>, T>;

template <class T>
inline constexpr bool is_view = detray::ranges::range<T>::value and
    std::is_object_v<T>&& std::is_move_constructible_v<T>and enable_view<T>;

/// Simply import the std versions of basic iterator functionality
/// @{
using std::begin;
using std::cbegin;
using std::cend;
using std::end;
using std::rbegin;
using std::rend;

using std::advance;
using std::distance;
using std::next;
using std::prev;
/// @}

/// Implementation of the detray ranges specific functions (size, data)
/// @{
using std::size;

/// @}

/// Compliance with std::ranges typedefs,
/// @see https://en.cppreference.com/w/cpp/ranges/iterator_t
/// @{
template <class T>
using iterator_t = decltype(detray::ranges::begin(std::declval<T&>()));

template <class T>
using const_iterator_t = iterator_t<const T>;

template <class T>
using sentinel_t = decltype(detray::ranges::end(std::declval<T&>()));

template <class T>
using range_size_t = decltype(detray::ranges::size(std::declval<T&>()));

template <class T>
using range_difference_t = typename std::iterator_traits<
    detray::ranges::iterator_t<T>>::difference_type;

template <class T>
using range_value_t =
    typename std::iterator_traits<detray::ranges::iterator_t<T>>::value_type;

template <class T>
using range_reference_t =
    typename std::iterator_traits<detray::ranges::iterator_t<T>>::reference;

template <class T>
using range_const_reference_t = const range_reference_t<T>;

template <class T>
using range_rvalue_reference_t = std::decay_t<range_value_t<T>>&&;
/// @}

// https://en.cppreference.com/w/cpp/ranges/view

/// Tags a type as a view
struct base_view {};

/// Defines a detray 'view'
template <typename view_impl_t>
struct view_interface : public base_view {

    constexpr view_interface() = default;

    DETRAY_HOST_DEVICE
    constexpr auto empty() const noexcept -> bool {
        return (m_impl_ptr->begin() == m_impl_ptr->end());
    }

    DETRAY_HOST_DEVICE
    constexpr auto data() noexcept { return &(*(m_impl_ptr->begin())); }

    DETRAY_HOST_DEVICE
    constexpr auto size() noexcept {
        return detray::ranges::distance(m_impl_ptr->begin(), m_impl_ptr->end());
    }

    DETRAY_HOST_DEVICE
    constexpr auto front() noexcept { return *(m_impl_ptr->begin()); }

    DETRAY_HOST_DEVICE
    constexpr auto back() noexcept {
        return *detray::ranges::next(m_impl_ptr->begin(), size() - 1);
    }

    DETRAY_HOST_DEVICE
    constexpr auto operator[](const dindex i) const {
        return *detray::ranges::next(m_impl_ptr->begin(), i);
    }

    view_impl_t* const m_impl_ptr{static_cast<view_impl_t*>(this)};
};

}  // namespace detray::ranges