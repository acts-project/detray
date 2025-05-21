/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <concepts>
#include <type_traits>

namespace detray::concepts {

template <typename T, typename... U>
concept any_of = std::disjunction_v<std::is_same<T, U>...>;

/// Arithmetic types
template <typename T>
concept arithmetic = std::is_arithmetic_v<T>;

template <typename T>
concept arithmetic_cvref = concepts::arithmetic<std::remove_cvref_t<T>>;

/// Same, except for cv qualifiers and reference
template <typename T, typename U>
concept same_as_cvref =
    std::same_as<std::remove_cvref_t<T>, std::remove_cvref_t<U>>;

/// Concept for detecting when a type is a non-const version of another
template <typename T, typename U>
concept same_as_no_const = std::same_as<std::remove_cv_t<T>, U>;

/// Index concept to access vector/matrix elements
template <typename T>
concept index = std::is_integral_v<T> && !std::same_as<T, bool>;

/// Concept that checks if a type models an interval of some value that can
/// be obtained with 'get'.
template <typename I>
concept interval = requires(I i) {

    requires(!concepts::arithmetic_cvref<I>);

    { detray::detail::get<0>(i) }
    ->concepts::arithmetic_cvref;

    { detray::detail::get<1>(i) }
    ->concepts::arithmetic_cvref;
};

}  // namespace detray::concepts
