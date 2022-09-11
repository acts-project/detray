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
#include "detray/utils/ranges/detail/iota_iterator.hpp"

namespace std {

/// Specializations of the std::iterator_traits struct for detray types
template <typename T>
struct iterator_traits<detray::ranges::detail::iota_iterator<T>> {
    using difference_type = T;
    using value_type = T;
    using pointer = T *;
    using reference = T &;
    using iterator_category = std::forward_iterator_tag;
};

}  // namespace std