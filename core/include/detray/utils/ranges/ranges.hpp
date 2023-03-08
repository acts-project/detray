/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges/detail/iterator_functions.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <cassert>
#include <memory>
#include <type_traits>

namespace detray::ranges {

/// @brief Provides c++17 detray iterators in a simplified std::ranges style,
///        meant to be used in device code.
///
/// @note Does make use of concepts and des not implement full ranges standard
/// (e.g. not every requirement is modelled, range/view complexity guarantees
/// are not strictly given and there is no lazy-evaluation).
///
/// @see https://en.cppreference.com/w/cpp/ranges
/// @see https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p0896r4.pdf
/// @{

/// Pull in iterator functionality
/// @{
// std version
using detray::ranges::detail::begin;
using detray::ranges::detail::cbegin;
using detray::ranges::detail::cend;
using detray::ranges::detail::crbegin;
using detray::ranges::detail::crend;
using detray::ranges::detail::end;
using detray::ranges::detail::rbegin;
using detray::ranges::detail::rend;

using detray::ranges::detail::data;
using detray::ranges::detail::empty;
using detray::ranges::detail::size;

// These are the only ones that are reimplemented using device quilifiers
using detray::ranges::detail::advance;
using detray::ranges::detail::distance;
using detray::ranges::detail::next;
using detray::ranges::detail::prev;
/// @}

/// Ranges typedefs,
/// @see https://en.cppreference.com/w/cpp/ranges/iterator_t
/// @{
template <class R>
using iterator_t = decltype(detray::ranges::begin(std::declval<R&>()));

template <class R>
using sentinel_t = decltype(detray::ranges::end(std::declval<R&>()));

template <class R>
using const_iterator_t = decltype(
    detray::ranges::begin(std::declval<const std::remove_reference_t<R>&>()));

template <class R>
using range_size_t = decltype(detray::ranges::size(std::declval<R&>()));

template <class R>
using range_difference_t =
    typename std::iterator_traits<detray::ranges::iterator_t<
        detray::detail::remove_cvref_t<R>>>::difference_type;

template <class R>
using range_value_t = typename std::iterator_traits<
    detray::ranges::iterator_t<detray::detail::remove_cvref_t<R>>>::value_type;

template <class R>
using range_reference_t = typename std::iterator_traits<
    detray::ranges::iterator_t<detray::detail::remove_cvref_t<R>>>::reference;

template <class R>
using range_const_reference_t = const range_reference_t<R>;

template <class R>
using range_rvalue_reference_t = std::add_rvalue_reference_t<range_value_t<R>>;

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
struct range<R,
             std::enable_if_t<
                 (std::is_class_v<detray::ranges::iterator_t<R>> or
                  std::is_pointer_v<detray::ranges::iterator_t<
                      R>>)and(std::is_class_v<detray::ranges::sentinel_t<R>> or
                              std::is_pointer_v<detray::ranges::sentinel_t<R>>),
                 void>> : public std::true_type {};

// Range
template <class R>
inline constexpr bool range_v = detray::ranges::range<R>::value;

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

// Range categories
template <class R>
inline constexpr bool input_range_v = detray::ranges::range_v<R>and
    input_iterator_v<detray::ranges::iterator_t<R>>;

template <class R>
inline constexpr bool output_range_v = detray::ranges::range_v<R>and
    output_iterator_v<detray::ranges::iterator_t<R>>;

template <class R>
inline constexpr bool forward_range_v = detray::ranges::range_v<R>and
    forward_iterator_v<detray::ranges::iterator_t<R>>;

template <class R>
inline constexpr bool bidirectional_range_v = detray::ranges::range_v<R>and
    bidirectional_iterator_v<detray::ranges::iterator_t<R>>;

template <class R>
inline constexpr bool random_access_range_v = detray::ranges::range_v<R>and
    random_access_iterator_v<detray::ranges::iterator_t<R>>;

// Contiguous iterator trait is only available in c++20
/*
template <typename R>
inline constexpr bool contiguous_range_v = ...
*/

/// @see https://en.cppreference.com/w/cpp/ranges/sized_range
template <class R>
inline constexpr bool disable_sized_range = false;

/// @brief A function 'size' is implemented for the range @tparam R
template <class R>
inline constexpr bool sized_range =
    not ranges::disable_sized_range<detray::detail::remove_cvref_t<R>> and
    (detray::ranges::range_v<R> and
     std::is_integral_v<detray::ranges::range_size_t<R>>);

/// @see https://en.cppreference.com/w/cpp/ranges/borrowed_range
template <class R>
inline constexpr bool enable_borrowed_range = false;

template <class R>
inline constexpr bool borrowed_range =
    detray::ranges::range_v<R> &&
    (std::is_lvalue_reference_v<R> ||
     ranges::enable_borrowed_range<detray::detail::remove_cvref_t<R>>);

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

/// Defines a detray 'view'. For now, the views have to restrict the member
/// functions themselves (i.e. forward, bidirectional, random access)
template <typename view_impl_t>
class view_interface : public base_view {

    public:
    constexpr view_interface() = default;

    /// @note requires forward range
    template <typename R = view_impl_t,
              std::enable_if_t<detray::ranges::forward_range_v<R>, bool> = true>
    DETRAY_HOST_DEVICE constexpr auto empty() const -> bool {
        return (_impl_ptr->begin() == _impl_ptr->end());
    }

    DETRAY_HOST_DEVICE
    constexpr explicit operator bool() const {
        return not detray::ranges::empty(*_impl_ptr);
    }

    /// @note requires contiguous range (not yet modelled)
    DETRAY_HOST_DEVICE
    constexpr auto data() const {
        return empty() ? nullptr : std::addressof(*(_impl_ptr->begin()));
    }

    /// @note requires contiguous range (not yet modelled)
    DETRAY_HOST_DEVICE
    constexpr auto data() {
        return empty() ? nullptr : &(*(_impl_ptr->begin()));
    }

    /// @note requires forward range
    template <typename R = view_impl_t,
              std::enable_if_t<detray::ranges::forward_range_v<R>, bool> = true>
    DETRAY_HOST_DEVICE constexpr auto size() const {
        return detray::ranges::distance(_impl_ptr->begin(), _impl_ptr->end());
    }

    /// @note requires forward range
    template <typename R = view_impl_t,
              std::enable_if_t<detray::ranges::forward_range_v<R>, bool> = true>
    DETRAY_HOST_DEVICE constexpr decltype(auto) front() const {
        const auto bg = _impl_ptr->begin();
        assert(not empty());
        return *bg;
    }

    /// @note requires forward range
    template <typename R = view_impl_t,
              std::enable_if_t<detray::ranges::forward_range_v<R>, bool> = true>
    DETRAY_HOST_DEVICE constexpr decltype(auto) front() {
        const auto bg = _impl_ptr->begin();
        assert(not empty());
        return *bg;
    }

    /// @note requires bidirectional range
    template <
        typename R = view_impl_t,
        std::enable_if_t<detray::ranges::bidirectional_range_v<R>, bool> = true>
    DETRAY_HOST_DEVICE constexpr decltype(auto) back() const {
        auto sentinel = _impl_ptr->end();
        assert(not empty());
        return *(--sentinel);
    }

    /// @note requires bidirectional range
    template <
        typename R = view_impl_t,
        std::enable_if_t<detray::ranges::bidirectional_range_v<R>, bool> = true>
    DETRAY_HOST_DEVICE constexpr decltype(auto) back() {
        auto sentinel = _impl_ptr->end();
        assert(not empty());
        return *(--sentinel);
    }

    /// Subscript operator that takes detray @c dindex
    ///
    /// @note requires random access range
    template <
        typename R = view_impl_t,
        std::enable_if_t<detray::ranges::random_access_range_v<R>, bool> = true>
    DETRAY_HOST_DEVICE constexpr decltype(auto) operator[](
        const dindex i) const {
        return (_impl_ptr->begin())
            [static_cast<detray::ranges::range_difference_t<R>>(i)];
    }
    /// Subscript operator that takes detray @c dindex
    ///
    /// @note requires random access range
    template <
        typename R = view_impl_t,
        std::enable_if_t<detray::ranges::random_access_range_v<R>, bool> = true>
    DETRAY_HOST_DEVICE constexpr decltype(auto) operator[](const dindex i) {
        return (_impl_ptr->begin())
            [static_cast<detray::ranges::range_difference_t<R>>(i)];
    }

    private:
    view_impl_t* const _impl_ptr{static_cast<view_impl_t*>(this)};
};

/// View traits
/// @{
template <typename R>
inline constexpr bool enable_view =
    std::is_base_of_v<base_view, R> or std::is_base_of_v<view_interface<R>, R>;

template <class R>
inline constexpr bool view = detray::ranges::range_v<R>and std::is_object_v<R>&&
    std::is_move_constructible_v<R>and enable_view<R>;

template <class R>
inline constexpr bool viewable_range =
    detray::ranges::range_v<R> &&
    (borrowed_range<R> || view<detray::detail::remove_cvref_t<R>>);
/// @}

}  // namespace detray::ranges