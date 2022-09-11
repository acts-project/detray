/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <type_traits>

#include "detray/definitions/detail/accessor.hpp"

namespace detray::detail {

/// Helper trait for detecting when a type is a non-const version of another
///
/// This comes into play multiple times to enable certain constructors
/// conditionally through SFINAE.
///
template <typename CTYPE, typename NCTYPE>
struct is_same_nc {
    static constexpr bool value = false;
};

template <typename TYPE>
struct is_same_nc<const TYPE, TYPE> {
    static constexpr bool value = true;
};

/// Helper trait that checks if a type models an interval of some value that can
/// be obtained with 'get'.
template <typename TYPE, typename = void, typename = void>
struct is_interval : public std::false_type {};

template <typename TYPE>
struct is_interval<
    TYPE,
    std::enable_if_t<
        not std::is_arithmetic_v<std::remove_reference_t<TYPE>> and
            std::is_arithmetic_v<std::remove_reference_t<
                decltype(detray::detail::get<0>(std::declval<TYPE>()))>>,
        void>,
    std::enable_if_t<
        not std::is_arithmetic_v<std::remove_reference_t<TYPE>> and
            std::is_arithmetic_v<std::remove_reference_t<
                decltype(detray::detail::get<1>(std::declval<TYPE>()))>>,
        void>> : public std::true_type {};

template <typename TYPE>
inline constexpr bool is_interval_v = is_interval<TYPE>::value;

}  // namespace detray::detail
