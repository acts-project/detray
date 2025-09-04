/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"

// System include(s)
#include <type_traits>

namespace detray::detail {

/// Helper trait that should emulate std::remove_cvref_t in c++20
template <typename T>
using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

/// Helper trait for detecting when a type is a non-const version of another
///
/// This comes into play multiple times to enable certain constructors
/// conditionally through SFINAE.
///
/// @{
template <typename CTYPE, typename NCTYPE>
struct is_same_nc {
    static constexpr bool value = false;
};

template <typename TYPE>
struct is_same_nc<const TYPE, TYPE> {
    static constexpr bool value = true;
};
/// @}

/// Extract the value type of e.g. a container.
/// @{
template <typename T, typename = void>
struct get_value_type {
    using type = T;
};

template <typename T>
struct get_value_type<T*, void> {
    using type = T;
};

template <typename container_t>
struct get_value_type<
    container_t,
    std::enable_if_t<
        not std::is_same<typename remove_cvref_t<container_t>::value_type,
                           void>::value ,
        void>> {
    using type = typename remove_cvref_t<container_t>::value_type;
};

template <typename T>
using get_value_t = typename get_value_type<T>::type;
/// @}

/// Helper trait that checks if a type models an interval of some value that can
/// be obtained with 'get'.
/// @{
template <typename TYPE, typename = void, typename = void>
struct is_interval : public std::false_type {};

template <typename TYPE>
struct is_interval<
    TYPE,
    std::enable_if_t<not std::is_arithmetic<std::remove_reference_t<TYPE>>::value  and
                         std::is_arithmetic<std::remove_reference_t<decltype(
                             detray::detail::get<0>(std::declval<TYPE>()))>>::value ,
                     void>,
    std::enable_if_t<not std::is_arithmetic<std::remove_reference_t<TYPE>>::value  and
                         std::is_arithmetic<std::remove_reference_t<decltype(
                             detray::detail::get<1>(std::declval<TYPE>()))>>::value ,
                     void>> : public std::true_type {};

template <typename TYPE>
inline constexpr bool is_interval_v = is_interval<TYPE>::value;
/// @}

/// Extract the first type from a parameter pack.
/// @{
template <typename first_t = void, typename... other_t>
struct first_type {
    using type = first_t;
};

template <typename... TYPES>
using first_t = typename first_type<TYPES...>::type;

/// Extract the first value from an index pack.
template <std::size_t first_t, std::size_t... other_t>
struct first_idx {
    static constexpr std::size_t value = first_t;
};

template <std::size_t... TYPES>
inline constexpr std::size_t first_idx_v = first_idx<TYPES...>::value;
/// @}

/// Get the position of a type in a template parameter pack
/// @{
template <typename T, typename... Ts>
struct get_type_pos {

    /// Unroll a parameter pack without using a tuple.
    ///
    /// @note Returns the position of the type counted from the back!
    template <typename first_t, typename... remaining_types>
    DETRAY_HOST_DEVICE static inline constexpr std::size_t type_pos_back() {
        if constexpr (not std::is_same<T, first_t>::value ) {
            return type_pos_back<remaining_types...>();
        }
        if constexpr (std::is_same<T, first_t>::value ) {
            return sizeof...(remaining_types) + 1;
        }
        return std::numeric_limits<std::size_t>::max();
    }

    static constexpr std::size_t value = sizeof...(Ts) - type_pos_back<Ts...>();
};

template <typename T, typename... Ts>
inline constexpr std::size_t get_type_pos_v = get_type_pos<T, Ts...>::value;
/// @}

template <class grid_t>
struct is_grid : public std::false_type {};

template <typename T>
inline constexpr bool is_grid_v = is_grid<T>::value;

template <class accelerator_t, typename = void>
struct is_surface_grid : public std::false_type {};

template <typename T>
inline constexpr bool is_surface_grid_v = is_surface_grid<T>::value;

template <class material_t, typename = void>
struct is_hom_material : public std::false_type {};

template <typename T>
inline constexpr bool is_hom_material_v = is_hom_material<T>::value;

template <class material_t, typename = void>
struct is_material_map : public std::false_type {};

template <typename T>
inline constexpr bool is_material_map_v = is_material_map<T>::value;

template <class material_t, typename = void>
struct is_volume_material : public std::false_type {};

template <typename T>
inline constexpr bool is_volume_material_v = is_volume_material<T>::value;

template <class material_t, typename = void>
struct is_surface_material : public std::false_type {};

template <typename T>
inline constexpr bool is_surface_material_v = is_surface_material<T>::value;

}  // namespace detray::detail
