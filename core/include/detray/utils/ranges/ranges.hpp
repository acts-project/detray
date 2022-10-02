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

/// @brief Provides c++17 detray iterators in std::ranges style, meant to be
////       used in device code.
///
/// @note Does make use of concepts and des not implement full ranges standard
/// (e.g. not every requirement is modelled, range/view complexity guarantees
/// are not given and there is no lazy-evaluation).
/// missing: ranges::sized_range, ranges::viewable_range, ranges::constant_range
///
/// @see https://en.cppreference.com/w/cpp/ranges
/// @{

/// Simply import the std versions of basic iterator functionality for now
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

using std::advance;
using std::distance;
using std::next;
using std::prev;
/// @}

/// Ranges typedefs,
/// @see https://en.cppreference.com/w/cpp/ranges/iterator_t
/// @{
template <class R>
using iterator_t = decltype(detray::ranges::begin(std::declval<R&>()));

template <class R>
using const_iterator_t = iterator_t<const R>;

template <class R>
using sentinel_t = decltype(detray::ranges::end(std::declval<R&>()));

template <class R>
using range_size_t = decltype(detray::ranges::size(std::declval<R&>()));

template <class R>
using range_difference_t = typename std::iterator_traits<
    detray::ranges::iterator_t<R>>::difference_type;

template <class R>
using range_value_t =
    typename std::iterator_traits<detray::ranges::iterator_t<R>>::value_type;

template <class R>
using range_reference_t =
    typename std::iterator_traits<detray::ranges::iterator_t<R>>::reference;

template <class R>
using range_const_reference_t = const range_reference_t<R>;

template <class R>
using range_rvalue_reference_t = std::decay_t<range_value_t<R>>&&;

/// @}

/// Range traits
/// @{

/// Base case: The type is not a 'range'
template <class R, typename = void>
struct range : public std::false_type {
    using type = void;
};

/// @brief A range has a begin() and an end() member function
///
/// @see https://en.cppreference.com/w/cpp/ranges/range
///
/// Existance of 'begin' and 'end' is guranteed by simply calling std::begin
/// and std::end.
/// @note In case of @c vecmem::device_vector the iterator is a pointer type.
template <class R>
struct range<
    R, std::enable_if_t<
           (std::is_class_v<std::decay_t<decltype(detray::ranges::begin(
                std::declval<R&>()))>> or
            std::is_pointer_v<std::decay_t<decltype(detray::ranges::begin(
                std::declval<R&>()))>>)and(std::
                                               is_class_v<std::decay_t<
                                                   decltype(detray::ranges::end(
                                                       std::declval<R&>()))>> or
                                           std::is_pointer_v<std::decay_t<
                                               decltype(detray::ranges::end(
                                                   std::declval<R&>()))>>),
           void>> : public std::true_type {};

template <class R>
inline constexpr bool range_v = detray::ranges::range<R>::value;

template <class R>
inline constexpr bool input_range_v = detray::ranges::range_v<R>and
    std::is_base_of_v<std::input_iterator_tag,
                      typename std::iterator_traits<
                          detray::ranges::iterator_t<R>>::iterator_category>;

template <class R>
inline constexpr bool output_range_v = detray::ranges::range_v<R>and
    std::is_base_of_v<std::output_iterator_tag,
                      typename std::iterator_traits<
                          detray::ranges::iterator_t<R>>::iterator_category>;

template <class R>
inline constexpr bool forward_range_v = detray::ranges::range_v<R>and
    std::is_base_of_v<std::forward_iterator_tag,
                      typename std::iterator_traits<
                          detray::ranges::iterator_t<R>>::iterator_category>;

template <class R>
inline constexpr bool bidirectional_range_v = detray::ranges::range_v<R>and
    std::is_base_of_v<std::bidirectional_iterator_tag,
                      typename std::iterator_traits<
                          detray::ranges::iterator_t<R>>::iterator_category>;

template <class R>
inline constexpr bool random_access_range_v = detray::ranges::range_v<R>and
    std::is_base_of_v<std::random_access_iterator_tag,
                      typename std::iterator_traits<
                          detray::ranges::iterator_t<R>>::iterator_category>;

// Available only in c++20
/*template <typename R>
inline constexpr bool contiguous_range_v = range_v<R> and
std::is_base_of_v<std::contiguous_iterator_tag, typename
std::iterator_traits<detray::ranges::iterator_t<R>>::iterator_category>;*/

/// @see https://en.cppreference.com/w/cpp/ranges/borrowed_range
template <class R>
inline constexpr bool enable_borrowed_range = false;

template <class R>
inline constexpr bool borrowed_range =
    detray::ranges::range_v<R> &&
    (std::is_lvalue_reference_v<R> ||
     ranges::enable_borrowed_range<
         std::remove_reference_t<std::remove_cv_t<R>>>);

/// @brief models a dangling iterator
/// @see https://en.cppreference.com/w/cpp/ranges/dangling
struct dangling {
    constexpr dangling() noexcept = default;
    template <class... Args>
    constexpr dangling(Args&&...) noexcept {}
};

template <class R>
using borrowed_iterator_t =
    std::conditional_t<borrowed_range<R>, detray::ranges::iterator_t<R>,
                       dangling>;

/// @see https://en.cppreference.com/w/cpp/ranges/common_range
template <class R>
inline constexpr bool common_range =
    detray::ranges::range_v<R>&& std::is_same_v<detray::ranges::iterator_t<R>,
                                                detray::ranges::sentinel_t<R>>;
/// @}

/// Definition of 'view'
/// @see https://en.cppreference.com/w/cpp/ranges/view
/// @{

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
        const auto i = size();
        return *detray::ranges::next(m_impl_ptr->begin(), i - decltype(i){1});
    }

    DETRAY_HOST_DEVICE
    constexpr auto operator[](const dindex i) const {
        return *detray::ranges::next(m_impl_ptr->begin(), i);
    }

    view_impl_t* const m_impl_ptr{static_cast<view_impl_t*>(this)};
};

/// View traits
/// @{
template <typename R>
inline constexpr bool enable_view =
    std::is_base_of_v<base_view, R> or std::is_base_of_v<view_interface<R>, R>;

template <class R>
inline constexpr bool view = detray::ranges::range_v<R>and std::is_object_v<R>&&
    std::is_move_constructible_v<R>and enable_view<R>;
/// @}

}  // namespace detray::ranges