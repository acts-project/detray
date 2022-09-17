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
#include "detray/utils/ranges/detail/chain_iterator.hpp"
#include "detray/utils/ranges/detail/enumerate_iterator.hpp"
#include "detray/utils/ranges/detail/iota_iterator.hpp"
#include "detray/utils/ranges/detail/single_iterator.hpp"

namespace std {

/// Specializations of std::iterator_traits struct for the iota range factory
template <typename T>
struct iterator_traits<detray::ranges::detail::iota_iterator<T>> {
    using difference_type = T;
    using value_type = T;
    using pointer = T *;
    using reference = T &;
    using iterator_category = std::input_iterator_tag;
};

/// Specializations of std::iterator_traits struct for the enumrator
template <typename T, typename I>
struct iterator_traits<detray::ranges::detail::enumerate_iterator<T, I>> {
    using difference_type = typename std::iterator_traits<T>::difference_type;
    using value_type = typename std::iterator_traits<T>::value_type;
    using pointer = typename std::iterator_traits<T>::pointer;
    using reference = typename std::iterator_traits<T>::reference;
    using iterator_category = std::input_iterator_tag;
};

/// Specializations of std::iterator_traits struct for the a single value
template <typename T>
struct iterator_traits<detray::ranges::detail::single_iterator<T>> {
    using difference_type = std::size_t;
    using value_type = T;
    using pointer = T *;
    using reference = T &;
    using iterator_category = std::bidirectional_iterator_tag;
};

/// Specializations of std::iterator_traits struct for the sequential iteration
/// of multiple ranges
template <std::size_t I, typename T>
struct iterator_traits<detray::ranges::detail::chain_iterator<I, T>> {
    private:
    using iterator_t = decltype(detray::detail::get<0>(std::declval<T>()));

    public:
    using difference_type =
        typename std::iterator_traits<iterator_t>::difference_type;
    using value_type = typename std::iterator_traits<iterator_t>::value_type;
    using pointer = typename std::iterator_traits<iterator_t>::pointer;
    using reference = typename std::iterator_traits<iterator_t>::reference;
    using iterator_category = std::forward_iterator_tag;
};

/// Specializations of std::iterator_traits struct for concurrent iteration of
/// multiple ranges
/*template <std::size_t I, typename T>
struct iterator_traits<detray::ranges::detail::zip_iterator<I, T>> {
    using difference_type = std::size_t;
    using value_type = T;
    using pointer = T *;
    using reference = T &;
    using iterator_category = std::bidirectional_iterator_tag;
};*/

}  // namespace std