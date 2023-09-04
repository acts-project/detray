/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"

// System include(s)
#include <type_traits>

namespace detray {

namespace detail {
/// The detray scalar/boolean types (can be SIMD)
/// @{
template <typename T, typename = void>
struct get_bool {};

template <typename T>
struct get_bool<T, std::enable_if_t<std::is_arithmetic_v<T>, void>> {
    using boolean = bool;
};

template <typename T>
struct get_bool<
    T, std::enable_if_t<!std::is_same_v<typename T::scalar, void>, void>> {
    using boolean = typename T::boolean;
};
/// @}

}  // namespace detail

template <typename A>
using dbool = typename detail::get_bool<A>::boolean;

namespace detail {

/// boolean utilities
/// @{
inline constexpr bool any_of(bool b) {
    // Do nothing
    return b;
}

#if (IS_SOA)
template <typename T,
          std::enable_if_t<Vc::Traits::is_simd_mask<T>::value, bool> = true>
inline bool any_of(T &&mask) {
    return Vc::any_of(std::forward<T>(mask));
}
#endif
/// @}

}  // namespace detail

}  // namespace detray
